#!/usr/bin/env python3

import string, random, optparse
from abc import ABC, abstractmethod

#This script reads in a fasta file and creates two outputs
    # 1) A new version of the fasta with names replaces by codes
    # 2) A text file linking the original names to the new coded names
# Additionally, this script can decode a supplied fasta (or fasta epitope map)

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option( '-f', '--fasta',  help='Fasta file. [None, REQ]' )
    p.add_option( '-l', '--length',  help='Length of code to use. If not provided, an optimal length will be calculated based on # of seqs [None]' )
    p.add_option( '-d', '--decode', help = "Use this flag if you wish to decode an already coded fasta ( or an epitope map ) "
                                           "containing coded names. Include with this flag the name of the coded file you want decoded, the file will be "
                                           "treated as either a coded epitope map, or a coded fasta file (this will be detected automatically). "
                                           "If this flag is included, then the '--key' flag must also be included."
                                           "When this flag is included, the decoded file will be written to the base name "
                                           "of '--key', with _decoded appended.",
                )
    p.add_option( '-k', '--key', help =  "Key mapping a coded name to its original name, this file is provided "
                                         "as output when coding a fasta file using this script. Note that this "
                                         "argument will only be considered when the '--decode' flag is also supplied."
                )

    opts, args = p.parse_args()
    
    if not opts.decode:
        encode_fasta( opts )
    else:
        decode_file( opts )

#----------------------End of main()

def write_code(names, codes, outfile):
    fout=open(outfile, "w")
    for i in range(len(names)):
        fout.write("%s\t%s\n" % (names[i], codes[i]))
    fout.close()

def rand_code(chars, l):
    code=""
    for each in range(l):
        code+=random.choice(chars)
#    print code
    return code

def read_fasta_dict_upper(file):
    names, seqs = read_fasta_lists(file)
    seqs = [x.upper() for x in seqs]
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                # Only grab the first item in the line,
                # protects against SEQ followed by anything else
                seqs.append( seq.strip().split()[ 0 ] )
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()

def encode_fasta( opts ):
        #Read in fasta file
        names, seqs = read_fasta_lists(opts.fasta)
        
        #Characters to use in code
        chars=string.ascii_letters + string.digits
        
        #Determine the size of code to use, if not provided
        if not opts.length:
            length=1
            while len(chars)**length<(2*len(names)):
                length+=1
        
        #Generate codes
        codes = set()
        while len(codes)<len(names):
            codes.add(rand_code(chars,length))
        codes=list(codes)
        
        write_fasta(codes, seqs, "coded_%s" % opts.fasta)
        write_code(names, codes, "coded_%s_key.txt" % (".".join(opts.fasta.split(".")[:-1])))

def decode_file( opts ):
    if is_fasta( opts.decode ):
        decoder = DecoderFactory().create_decoder( "fasta_decoder" )
    else:
        decoder = DecoderFactory().create_decoder( "map_decoder" )

    decoder.set_key_file( opts.key )
    decoder.read_key()
    decoder.decode( opts.decode )
    decoder.write_output( "%s_decoded" % opts.decode )

def is_fasta( filename ):
    looking_for_name = True
    with open( filename, 'r' ) as open_file:
        line = open_file.readline()
        # In a fasta file, the first line should start with '>'
        if line.strip()[ 0 ] != '>':
            return False
        return True
    
class FileDecoder( ABC ):

    def __init__( self, filename = None ):
        self._key_file     = filename
        self._key          = {}

        self._decoded_data = None
        self._data         = None

    @abstractmethod
    def decode( self, filename ):
        pass

    @abstractmethod
    def write_output( self, filename ):
        pass

    def set_key_file( self, filename ):
        self._key_file = filename

    def read_key( self ):
        if not self._key:
            self._key = self._parse_key( self._key_file )
        return self._key

    def _parse_key( self, filename ):
        parsed_dict = {}
        with open( filename, 'r' ) as open_key:
            for line in open_key:
                split_line = line.split( '\t' )
                coded_name   = split_line[ 1 ].strip()
                uncoded_name = split_line[ 0 ].strip()

                parsed_dict[ coded_name ] = uncoded_name
        return parsed_dict

class FastaDecoder( FileDecoder ):
    def __init__( self, filename = None ):
        super().__init__( filename = filename )

    def decode( self, filename ):
        coded_fasta  = self._fasta_to_dict( filename )
        key_dict     = self.read_key()
        decoded_data = {}

        for coded_name, sequence in coded_fasta.items():
            decoded_data[ key_dict[ coded_name ] ] = sequence

        self._decoded_data = decoded_data

    def write_output( self, filename ):
        self._write_fasta( filename )

    def _write_fasta( self, outfile_name ):
        with open( outfile_name, 'w' ) as out_file:
            for name, sequence in self._decoded_data.items():
                out_file.write( ">%s\n%s\n" % ( name, sequence ) )

        
    def _fasta_to_dict( self, filename ):
        out_dict = {}
        current_seq = [ "", "" ] # List containing ( name, sequence )
        with open( filename, 'r' ) as fasta:
            for line in fasta:
                if line[ 0 ] == '>':
                    # Note that dict will contain an entry { "":""}
                    out_dict[ current_seq[ 0 ] ] = current_seq[ 1 ] 
                    current_seq = [ line.strip().strip( ">" ), "" ]
                else:
                    current_seq[ 1 ] += line.strip()

        # Delete the { "": "" } entry
        del out_dict[ '' ]

        return out_dict

class MapDecoder( FileDecoder ):
    def __init__( self, filename = None ):
        super().__init__( filename = filename )
        self._delimiter    = '~'
        self._decoded_data = {}
        self._open_decode_file = None

    def decode( self, filename ):
        self.read_key()

        if self._open_decode_file is None:
            self._open_decode_file = open( filename, 'r' )
                        
    def _decode_line( self, line_data ):
        DELIMITER = self._delimiter
        out_data  = {}
        split_data = line_data.split( '\t' )

        key               = split_data[ 0 ].strip()
        decoded_data_list = list()

        split_data = split_data[ 1 ].split( DELIMITER )

        for item in split_data:
            item = self._decode_item( item.strip() )
            decoded_data_list.append( item )
        return [ key, decoded_data_list ]

    def _get_decoded_data( self ):
        # Decode more data, one line from the file at a time
        while True:
            current_data = self._open_decode_file.readline()
            if current_data:
                decoded_line = self._decode_line( current_data )
                yield decoded_line
            else:
                break
                    

    def _decode_item( self, item ):
        return self._key[ item ]

    def write_output( self, filename ):
        with open( filename, 'w' ) as open_file:
            for item in self._get_decoded_data():
                open_file.write( "%s\t%s\n" % ( item[ 0 ], self._delimiter.join( item[ 1 ] ) ) )
            

class DecoderFactory:
    MAP_DECODER   = "map_decoder"
    FASTA_DECODER = "fasta_decoder"

    def create_decoder( self, decoder_type, filename = None ):
        if decoder_type == DecoderFactory.MAP_DECODER:
            return MapDecoder( filename )
        elif decoder_type == DecoderFactory.FASTA_DECODER:
            return FastaDecoder( filename )

###------------------------------------->>>>    
if __name__ == "__main__":
    main()

