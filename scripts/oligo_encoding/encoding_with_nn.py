#!/usr/bin/env python3
import h2o      # For running models
import argparse # for parsing command-line args
import sys

def main():
    arg_parser = argparse.ArgumentParser( description = "Use h2o to select encodings for oligos." )  

    arg_parser.add_argument( '-m', '--model', description = "Trained model that ill be used to 

if __name__ == '__main__':
    main()
