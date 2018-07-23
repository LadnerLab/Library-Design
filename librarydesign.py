#!/usr/bin/env python3
import optparse
import sys
import time
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

    cluster_options = ( "-q %s -l %s -n %d -s %s -o %s -c %s --id %d -k %d"
                        % ( options.query, options.lineage, options.number, options.start, options.cluster_dir, options.cluster_method,
                            options.id, options.xmer_window_size
                          )  
                      ) 
    kmer_options = ( '-i %d  -x %d  -y %d -r %d -t %d' 
                     % ( options.iterations, options.xmer_window_size,
                         options.ymer_window_size, options.redundancy, options.threads
                       )
                   )
    if options.functional_groups:
        kmer_options += ' -p '
    if options.min_xmer_coverage:
        kmer_options += '-c ' + str( options.min_xmer_coverage )

    cluster_options += " -q combined.fasta "
    cluster_script = SBatchScript( "clustering.py " + cluster_options, "slurm_script",
                                   options.slurm
                                 )  

    cluster_script.add_module( "python/3.latest" )
    cluster_script.write_script()
    cluster_script.run()

    while not os.path.exists( options.cluster_dir ): 
        time.sleep( 1 )

    cluster_files = os.listdir( options.cluster_dir )
    num_files = len( cluster_files )

    while num_files < 1:
        time.sleep( 1 )
        num_files = len( os.listdir( options.cluster_dir ) )
    
    os.chdir( options.cluster_dir )
    job_ids = list()

    for current_file in cluster_files:
        kmer_options[ '-q' ] = current_file
        kmer_options[ '-o' ] = current_file + " out"

        kmer_script = SBatchScript( "kmer_oligo " + kmer_options, "kmer_script", options.slurm )
        kmer_script.write_script()
        output = str(kmer_script.run() ).split()
        job_ids.append( output[ -1 ] )

    job_ids = [ item.split( '\\n' )[ 0 ] for item in job_ids ]

    ids_combined = ",".join( job_ids )

    # combination_script_command = { "*_R_1 > " "combined.fasta; mv combined.fasta ../combined.fasta"  }
    # combination_script = SBatchScript( "cat ", "combine_script", combination_script_command, [ "dependency  afterany:" + ids_combined ]  )
    combination_script.write_script()
    combination_script.run()



        

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

class SBatchScript:
    def __init__( self, command, script, slurm_args, dependency_mode = "afterany" ):
        self.commands = [ SBatchScript.Command( command ) ]
        self.slurm_args = [ item.split() for item in slurm_args ]
        self.script = script

        self.dependencies = list()
        self.dependency_mode = dependency_mode

        self.modules = list()
        self.job_num = 0

        self.sbatch = "#SBATCH "
        self.shebang = "#!/bin/sh "

    class Command:
        def __init__( self, string_command ):
            self.command = string_command
        def __str__( self ):
            return self.command
        def add_arg( self, to_add ):
            self.command += to_add

    def write_script( self ):
        file = open( self.script, 'w' )

        file.write( self.shebang )
        file.write( "\n" )

        for item in self.slurm_args:
            file.write( self.sbatch + "--" + item[ 0 ] + "=" + item[ 1 ] )
            file.write( "\n" )

        if len( self.dependencies ) > 0:
            file.write( self.sbatch + "--dependency=" + self.dependency_mode + ','.join( self.dependencies ) )

        for current_module in self.modules:
            file.write( "module load " + current_module )
            file.write( "\n" )

        for current_command in self.commands:
            file.write( "srun " + str( current_command ) )
            file.write( "\n" )

        file.close()

    def run( self ):
        os.chmod( self.script, 0o755 )
        script = subprocess.getoutput( "sbatch " + self.script ) 

        # Get and return the jobnumber
        script = script.split()[ 3 ]
        self.job_num = script

        return script

    def add_command( self, in_command ):
        self.commands.append( SBatchScript.command( in_command ) )

    def add_slurm_arg( self, new_args ):
        for item in new_args:
            self.slurm_args.append( [ item.split() ] )

    def add_dependencies( self, job_num_list ):
        for current_job in nob_num_list:
            self.dependencies.append( current_job )

    def set_dependency_mode( self, new_mode ):
        self.dependency_mode = new_mode

    def add_modules( self, modules_list ):
        for item in modules_list:
            self.modules.append( item )
        
    def add_module( self, to_add ):
        self.modules.append( to_add )
       
        
            

  
if __name__ == '__main__':
    main()

