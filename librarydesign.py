#!/usr/bin/env python3
import optparse
import sys
import subprocess
import os

def main():
    usage = "usage %prog [options]"

    option_parser = optparse.OptionParser( usage ) 

    add_program_options( option_parser )

    options, arguments = option_parser.parse_args()

    check_required_option( options.query, "Fasta query file must be provided", True )
    if 'tax' in options.cluster_method:
        check_required_option( options.lineage, "Lineage file must be provided when using taxonomic clustering", True )

    cluster_options = { "-q": options.query, "-l": options.lineage, "-n": options.number, "-s": options.start,
                        "-o": options.cluster_dir, "-c": options.cluster_method, "--id": options.id 
                      } 
    kmer_options = { '-i': options.iterations, '-x': options.xmer_window_size, '-y': options.ymer_window_size,
                     '-r': options.redundancy, '-t': options.threads
                   }
    if options.functional_groups:
        kmer_options[ '-p' ] = ''
    if options.min_xmer_coverage:
        kmer_options[ '-c' ] = options.min_xmer_coverage

    cluster_script = SBatchScript( "clustering.py", "slurm_script", cluster_options,
                                   options.slurm
                                 )  

    cluster_script.write_script()
    cluster_script.run()

    cluster_files = os.listdir( options.cluster_dir )
    os.chdir( options.cluster_dir )

    for current_file in cluster_files:
        kmer_options[ '-q' ] = current_file
        kmer_options[ '-o' ] = current_file + " out"

        kmer_script = SBatchScript( "kmer_oligo", "kmer_script", kmer_options, options.slurm )
        kmer_script.write_script()
        kmer_script.run()


        

def add_program_options( option_parser ):
    option_parser.add_option( '-q', '--query', help = "Fasta query file to read sequences from and do ordering of. [None, Required]" )

    option_parser.add_option( '-l', '--lineage', help = "Taxonomic lineage file such as the one from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" )

    option_parser.add_option( '-n', '--number', type = int, default = 10000,
                              help = "Threshold value for determining cutoff of number of sequences that can be included in each output. [10,000]"
                            )
    option_parser.add_option( '-s', '--start', action = "append", default = 'family',
                              help = ( "Level of the taxonomic hierarchy at which to begin "
                                       "clustering. If this option is given multiple times, "
                                       "e.g. -s family -s phylum, "
                                       "they will be processed in order of taxonomic rank, e.g., "
                                       "superkingdom, kingdom, phylum, class, order, family, genus, species [ family ]"

                                     )
                            )
    option_parser.add_option( '-o', '--output', default = 'library.fasta',
                              help = "Name of oligo library file that will contain the final library design"
                            )
    option_parser.add_option( '-m', '--cluster_method', default = 'kmer',
                              help = ( "Method to use for clustering. Can be taxonomic or kmer-based. If taxonomic is selected, "
                                       "a taxonomic lineage file must also be provided. No lineage file is necessary for kmer "
                                       "clustering method. [kmer]"
                                     )
                            )
    option_parser.add_option( '--id', default = 0.8, type = float,
                              help = ( "Percentage of its kmers a sequence must share with a "
                                       "cluster in order for it to become a member of that cluster"
                                       "only used for kmer-based clustering [0.8]"
                                     )
                            )
    option_parser.add_option( '-x', '--xmer_window_size', type = 'int',
                              default = 10,
                              help = "Amount of characters from each Xmer alignment sequence to look at. [19]"
    )

    option_parser.add_option( '-y', '--ymer_window_size', type = 'int',
                              default = 19,
                              help = "Amount of characters from each Ymer alignment sequence to look at. [19]"
    )

    option_parser.add_option( '-r', '--redundancy', type = 'int', default = 1, help = "A number specifying the redundancy to be used to each kmer [1]" )

    option_parser.add_option( '-i', '--iterations', type = 'int', default = 1,
                              help = "Number of independent iterations to run. The result with the fewest oligos will be output [1]"
                            )

    option_parser.add_option( '-t', '--threads', type = 'int', default = 1,
                              help = "Number of threads to use when performing opterations [1]"
                            )

    option_parser.add_option( '-f', '--functional_groups', action = "store_true", dest = "functional_groups", default = False,
                              help = "Option to enable functional grouping of proteins"
                            )

    option_parser.add_option( '-c', '--min_xmer_coverage',
                              help = "Option to set the floating point minimum amount of coverage necessary for the program to cease execution. [1.0]"
                            )  

    option_parser.add_option( '--time', type = str,
                              help = "Time for given to each slurm script to run. Format is in days-hours:minutes:seconds, as specified by slurm. [1:00:00]"
                            )  

    option_parser.add_option( '--slurm', action = "append", 
                              help = ( 'slurm arguments to be written to the script, each should be entered as a separate ' 
                                       'argument such as: --slurm "mem 20G" --slurm "time 20:00" to specify a run with '
                                       '20 GB of memory that has 20 minutes '
                                     )  
                            )

    option_parser.add_option( '--cluster_dir', default = "tax_out",
                              help = "Name of directory to write clusters to. Note: this directory is created if it does not already exist. [tax_out]"
                            )



def check_required_option( option, string, exit_on_failure = False ):
    """
        Checks to see if a required option exists, prints out string and exits if that is not the case
    """
    if option is None:
        print( string )
        if exit_on_failure:
            print( "Exiting program due to above failures" )
            sys.exit( 0 )

def run_command_from_options( command_name, options_dict ):
    command = command_name + " "
    for flag, value in options_dict.items():
        command += str( flag )
        command += " "
        command += str( value )
        command += " "

    # Check that the script is in our path, or in the local directory
    script_found = script_exists( command_name )
    if script_found == "local":
        command = "./" + command
    elif not script_found:
        return False
    
    command = subprocess.Popen( "./" + command, shell = True )
    command.wait()

    return True
        
    
def script_exists( command_name ):
     file_found = True
     try:
         in_path = subprocess.check_output( [ command_name ] )
         file_found = "path"
     except FileNotFoundError:
         try:
             in_path = subprocess.check_output( [ "./" + command_name ] )
             file_found = "local"
         except FileNotFoundError:
             file_found = False
 
     return file_found

class SBatchScript():
    
    def __init__( self, command, output, program_args, slurm_args ):
        self.command = command 
        self.slurm_args = [ item.split() for item in slurm_args ]
        self.output = output
        self.program_args = program_args

        self.sbatch = "#SBATCH "
        self.shebang = "#!/bin/sh "

    def write_script( self ):
        file = open( self.output, 'w' )

        file.write( self.shebang )
        file.write( "\n" )

        for item in self.slurm_args:
            file.write( self.sbatch + "--" + item[ 0 ] + "=" + item[ 1 ] )
            file.write( "\n" )

        file.write( "module load python/3.latest\n" )
        file.write( "srun " + self.command + " " )

        for flag, argument in self.program_args.items():
            file.write( flag + " " + str( argument ) + " " )
        file.close()

    def run( self ):
        os.chmod( self.output, 0o755 )
        output = subprocess.Popen( "sbatch " + self.output, shell = True ) 
        output.wait()

    def add_program_arg( self, flag, arg ):
        self.program_program_args[ flag ] = arg
        
        
            

  
if __name__ == '__main__':
    main()

