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
    kmer_options = ( '-i %d  -x %d  -y %d -r %d -t %d ' 
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
        kmer_options += ' -q ' + str( current_file )
        kmer_options += ' -o ' + str( current_file ) + "_out "

        kmer_script = SBatchScript( "kmer_oligo " + kmer_options, "kmer_script", options.slurm )
        kmer_script.write_script()
        current_job_id = kmer_script.run()
        job_ids.append( current_job_id )

    combination_script = SBatchScript( "cat *_R_1 > combined.fasta", "combine_script", options.slurm, dependency_mode = "afterok" )
    combination_script.add_command( "mv combined.fasta ../combined.fasta" )
    combination_script.add_dependencies( job_ids )
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
    """
        Encapsulates a bash script to run on a server managed by slurm.
        Handles the writing and execution, of these scripts, and captures and stores
        any job numbers generated. This class also supports the import of modules,
        if this package is available on your system.
    """
    def __init__( self, command, script_name, slurm_args, dependency_mode = "afterany" ):
        """
            Constructor for SBatchScript class

            :param command: command to be run by slurm server, e.g., 'cat *.fasta'
             Note: bash shebang written to the file, but can be set by user if the standard
                  '#!/bin/sh' is not used by your system
             Note: srun will be prepended to the command, so the above becomes 'srun cat *.fasta'
             Note: multiple job steps can be included in a fasta file, but only one can be provided upon initialization
        
            :param script_name: name of the executable to be created by SBatchScript.write()

            :param slurm_args: list of slurm arguments to be written to the file. This
                               param is in the form of [ '--mem 4g', '--time 20:00', ... ]
             Note: the #SBATCH flag is written to the file before each of these arguments
        
            :param dependency_mode: Optional mode of dependencies this script is dependant upon.
        """
        self.commands = [ SBatchScript.Command( command ) ]
        self.slurm_args = [ item.split() for item in slurm_args ]
        self.script_name = script_name

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

    def write_script_name( self ):
        """
            Writes the script, the name of the executable created is 
            determined by the class-member variable script_name
        
            Note: this method sets the mode to octal 755 r/w access
        """
        file = open( self.script_name, 'w' )

        file.write( self.shebang )
        file.write( "\n" )

        for item in self.slurm_args:
            file.write( self.sbatch + "--" + item[ 0 ] + "=" + item[ 1 ] )
            file.write( "\n" )

        if len( self.dependencies ) > 0:
            file.write( self.sbatch + "--dependency=" + self.dependency_mode + ':' + ','.join( self.dependencies ) )
            file.write( "\n" )


        for current_module in self.modules:
            file.write( "module load " + current_module )
            file.write( "\n" )

        for current_command in self.commands:
            file.write( "srun " + str( current_command ) )
            file.write( "\n" )

        file.close()

    def run( self ):
        """
            Executes the script, and returns the slurm job number

            Note: this method sets the mode access mode to octal 755 
        """
        os.chmod( self.script_name, 0o755 )
        script = subprocess.getoutput( "sbatch " + self.script_name ) 

        # Get and return the jobnumber
        script = script.split()[ 3 ]
        self.job_num = script

        return script

    def set_shebang( self, new_shebang ):
        """
            Sets the shebang (Default '#!/bin/sh')
            to string new_shebang
        """
        self.shebang = new_shebang

    def add_command( self, in_command ):
        """
            Adds a command, (job-step) to be written to the output
            bash file.
        
            :param in_command: command to be written to the file, can be any
                               command recognized by your bash/slurm environment
            Note: srun is prepended to the command as it is written to the file, do not
                  include this yourself
        """
        self.commands.append( SBatchScript.Command( in_command ) )

    def add_slurm_arg( self, new_arg ):
        """
            Adds a new argument to be written to the executable created by this script,
        
            Note: before any slurm arguments are written to the file,
                  #SBATCH is written before any arguments, do not include it
                  here
            :param new_arg: argument to be written to script produced by this
                            obect's write method, in the format '--key value', or 
                            of the form '-c 1'
        """
        self.slurm_args.append( [ item.split() ] )

    def add_dependencies( self, job_num_list ):
        """
            Add a list of dependencies that this object relies upon.
            
            :param job_num_list: list of job numbers this script is to rely upon
        """
        for current_job in job_num_list:
            self.dependencies.append( current_job )

    def set_dependency_mode( self, new_mode ):
        """
            Sets the dependency mode of this job's dependencies
            :param new_mode: slurm dependent mode of dependencies,
                             can include 'afterany', 'afterok', etc
        """
        self.dependency_mode = new_mode

    def add_modules( self, modules_list ):
        """
            Adds a list of string modules to be loaded before execution of
            any job steps in the script. 
        
            Note: 'module load ' is written to the file by the script, do not include this
                  before any of the dependencies in modules_list
            :param modules_list: list of string modules to load
                                 [ 'python/3.6', 'blast+', ... ]
        """
        for item in modules_list:
            self.modules.append( item )
        
    def add_module( self, to_add ):
        """
            Add a single string module to the list of modules 
            that will be loaded before execution of any jobsteps in script.

            Note: 'module load ' is written to the file by the script, do not include this
                  in to_add 
        
            :param to_add: string module to add
        """
        self.modules.append( to_add )
       
        
            

  
if __name__ == '__main__':
    main()

