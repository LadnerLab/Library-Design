#!/usr/bin/env python3

import sys
import os
import argparse
import datetime

import matplotlib.pyplot as plt


def main():

    arg_parser = argparse.ArgumentParser( description = "Script to parse data produced by jobstats" )

    arg_parser.add_argument( '-s', '--stats', help = "File containing parseable jobstats output" )
    arg_parser.add_argument( '-k', '--kmers', help = "File containing kmer counts for each file run" )
    arg_parser.add_argument( '--num_yticks', default = 20, type = int )
    arg_parser.add_argument( '--suffix', help = "Suffix to remove from each name" )
    arg_parser.add_argument( '--max_yval', type = int, default = 9000 ) # multiple of 60
    arg_parser.add_argument( '--yaxis', help = "Data to show on yaxis", type = str, default = "elapsed" )
    arg_parser.add_argument( '--xaxis', help = "Data to show on xaxis", type = str, default = "kmers" )

    args = arg_parser.parse_args()

    file_parser = DataParser( args.stats )
    file_parser.read()

    runtime_data_from_file = file_parser.get_data()
    SUFFIX = args.suffix

    job_names = {}

    step_size = calc_step_size( args.max_yval, args.num_yticks )

    sizes_file = open( args.kmers, "r" )

    for line in sizes_file:
        line = line.split( '|' )
        job_names[ line[ 0 ] ] = line[ 1 ]

    job_array = list()
    if runtime_data_from_file:

        for current_data in runtime_data_from_file:
            data_job = SlurmJob( current_data )

            if data_job.is_completed() and SUFFIX in data_job._name:
                data_job.remove_suffix_from_name( SUFFIX )
                if data_job._name in job_names:
                    
                    data_job._kmers = job_names[ data_job._name ]
                    job_array.append( data_job )
            elif not data_job.is_completed() and SUFFIX in data_job._name:
                print( "Job %s (%s) completed with a data_jobstats of %s and was excluded" % ( data_job._name, data_job._id, data_job._state ) )
    else:
        print( "Unable to find data from statistics file, exiting..." )
        sys.exit( 1 )
                

    y_axis = get_axis_data( args.yaxis, job_array )
    x_axis = get_axis_data( args.xaxis, job_array )

    y_axis_label = get_axis_label( args.yaxis )
    x_axis_label = get_axis_label( args.xaxis )

    if y_axis is None or x_axis is None:
        print( "Incorrect specification for either the y or x axis, program cannot continue" )
        sys.exit( 1 )
            
    if len( y_axis ) > 0:
        y_tick_vals = get_yvals( y_axis, step_size )

        ax = plt.subplot()
        ax.set_yticks( y_tick_vals )
        ax.set_yticklabels( [ from_seconds( item ) for item in y_tick_vals ] )
        ax.plot()
        ax.scatter( x_axis, y_axis )
        plt.xlabel( "Number of kmers" )
        plt.ylabel( "Time (in seconds) to run " )
        plt.show()
    else:
        print( "No valid values were found to place on the y-axis" )

def to_seconds( string_time ):
    string_time = string_time.split( ':' )
    total = 0
    hours = int( string_time[ 0 ] )
    minutes = int( string_time[ 1 ] )
    seconds = int( string_time[ 2 ] )


    return ( 60 * hours * 60 ) + ( 60 * minutes ) + ( seconds )

def get_yvals( y_axis_vals, step_size ):
    return range( 0, max( y_axis_vals ) + 1, step_size )

def from_seconds( int_seconds ):
    hours = int_seconds // ( 60 * 60 )
    int_seconds %= ( 60 * 60 )
    minutes = int_seconds // 60
    int_seconds %= 60
    seconds = int_seconds

    date_time_obj = datetime.time( hour = hours, minute = minutes, second = seconds )
    return str( date_time_obj )

def to_datetime( string_time ):
    string_time = string_time.split( ':' )
    hours = int( string_time[ 0 ] )
    minutes = int( string_time[ 1 ] )
    seconds = int( string_time[ 2 ] )

    return datetime.time( hour = hours, minute = minutes, second = seconds )

def calc_step_size( max_val, num_ticks ):
    return  max_val // num_ticks

def get_axis_data( data_label, job_array ):
    return_array = None

    if data_label == "kmers":
        return_array = [ int( job._kmers ) for job in job_array ]
    elif data_label == "elapsed":
        return_array = [ to_seconds( job._elapsed ) for job in job_array ]
    elif data_label == "mem_used":
        return_array = [ int( job._used_mem ) for job in job_array ]
    elif data_label == "req_cpu":
        return_array = [ int( job._req_cpu ) for job in job_array ]
    elif data_label == "state":
        return_array = [ job._state for job in job_array ]
    elif data_label == "time_limit":
        return_array = [ to_seconds( job._time_limit ) for job in job_array ]
    elif data_label == "mem_req":
        return_array = [ job._req_mem for job in job_array ]

    return return_array

def get_axis_label( data_label ):
    return_str = None

    if data_label == "kmers":
        return_str = "Number of kmers"
    elif data_label == "elapsed":
        return_str = "Elapsed time (in minutes)"
    elif data_label == "mem_used":
        return_str = "Amount of Memory Used"
    elif data_label == "req_cpu":
        return_str = "Number of required CPUs"
    elif data_label == "state":
        return_str = "Job state"
    elif data_label == "time_limit":
        return_str = "Time limit (in minutes)"
    elif data_label == "mem_req":
        return_str = "Memory allocated for Job"

    return return_str
    


class DataParser:
    def __init__( self, input_file, delimiter_char = '|' ):

        self._input_file = input_file
        self._data = list()
        self._failed_to_read = False
        self._delimiter_char = delimiter_char

    def read( self ):
        try:
            read_file = open( self._input_file, "r" )

            for line in read_file:
                if line[ 0 ] != '#':
                    self._data.append( line.split( self._delimiter_char ) )

            read_file.close()

        except FileNotFoundError:
            print( "FILE %s NOT FOUND" % self._input_file )
            self._failed_to_read = True
    def get_data( self ):
        return self._data


class SlurmJob:
    def __init__( self, job_data_list ):
        self._job_data = job_data_list

        self._id = job_data_list[ 0 ]
        self._name = job_data_list[ 1 ]
        self._req_mem = job_data_list[ 2 ]

        if len( job_data_list[ 3 ] ) > 1:
            self._used_mem = job_data_list[ 3 ][:-1:]
        else:
            self._used_mem = job_data_list[ 3 ]
            
        self._req_cpu = job_data_list[ 4 ]
        self._user_cpu = job_data_list[ 5 ]
        self._time_limit = job_data_list[ 6 ]
        self._elapsed = job_data_list[ 7 ]
        self._state = job_data_list[ 8 ]
        self._kmers = ""

    def remove_suffix_from_name( self, suffix_to_remove ):
        new_name = self._name.split( suffix_to_remove )[ 0 ]
        self._name = new_name
        self._job_data[ 0 ] = new_name
    def is_completed( self ):
        return self._state == "COMPLETED"
       
if __name__ == '__main__':
    main()
