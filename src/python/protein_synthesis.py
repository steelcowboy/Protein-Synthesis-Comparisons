#!/usr/bin/env python
from synthesis import *
import argparse
from sys import exit

def main():
    parser = argparse.ArgumentParser(description='Translate DNA and compare against other RNA strands')
<<<<<<< HEAD
    parser.add_argument('-v', action='version', version=str(VERSION), 
=======
    parser.add_argument('-v', action='version', version='0.5', 
>>>>>>> origin/master
                        help='print the version')
    parser.add_argument('-o', dest='filename',
                        help='speficy a filename to output to')
    parser.add_argument('-t', action='store_true', dest='test_mode',
                        help='enable test mode which will run the program 10000 times')
    parser.add_argument('-q', action='store_true', dest='quiet',
                        help='do not print to stdout')
    parser.add_argument('-d', action='store_true', dest='default_mode', 
                        help='use a hardcoded example, mainly for testing')
    parser.add_argument('-r', action='store_true', dest='given_rna',
                        help='specify that RNA was given as input')
    parser.add_argument('dna', nargs=argparse.REMAINDER,
                        help='DNA strands to be processed')
    args = parser.parse_args()

    if len(argv) == 1:
        parser.print_help()
        exit(1)

    output_to_file = True if args.filename else False
    
    if args.default_mode:
        default_process(args.test_mode, args.quiet, output_to_file, args.filename)
    else:
        orgs = []
        numorgs = len(args.dna)
        
        if numorgs == 1:
            orgs.append(Organism("idk", 0, args.dna[0]))
        else:
            for i, x in enumerate(args.dna):
                name = input("What is strand {}? ".format(i+1))
                orgs.append(Organism(name, i, args.dna[i]))
            print("-"*20)
        outinfo = calculate_orgs(orgs, args.given_rna)
        output_info(outinfo, args.quiet, output_to_file, args.filename)


if __name__ == "__main__":
    main()
            
