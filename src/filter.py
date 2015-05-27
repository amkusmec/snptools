# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:47:00 2015

@author: aaron
"""

import sys
import argparse
import textwrap
import timeit
import os

###############################################################################
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

###############################################################################
def version():
   v0 = """
   ############################################################################
   filter
   (c) 2015 Aaron Kusmec
   
   Filter SNPs based on missing rates/minor allele frequencies.
   Input modes,
       1 = .dsf
       2 = .hmp.txt
       3 = .ped (PLINK)
       
   Usage: python3 filter.py -s example.stat -i example.dsf -o filtered -m 1 -n 0.6 -f 0.05
       
   ############################################################################
   """
   
   return v0

###############################################################################  
def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description = textwrap.dedent(version()))

    parser.add_argument('-p', '--path', help = 'Path of the input file', \
                        nargs = '?', default = os.getcwd())
    parser.add_argument('-s', '--stat', help = 'Stat file', type = str)
    parser.add_argument('-i', '--input', help = 'Input file', type = str)
    parser.add_argument('-o', '--output', help = 'Output file (no ext)', type = str)
    parser.add_argument('-m', '--mode', help = 'I/O mode', type = int)
    parser.add_argument('-n', '--miss', help = 'Max missing rate', \
                        type = float, default = 1.0)
    parser.add_argument('-f', '--maf', help = 'Minimum minor allele frequency',\
                        type = float, default = 0.0)

    return parser

###############################################################################
def checkFile(filename, mode):
    print("Checking [ ", filename, " ].")
    
    with open(filename, 'r') as infile:
        line = infile.readline().split()
    
    if mode == 1:
        if line[0] != "snpid" or line[1] != "major" or line[2] != "minor":
            warning(".dsf formatted incorrectly.")
    elif mode == 2:
        if line[0] != "rs" and line[0] != "rs#":
            warning(".hmp.txt formatted incorrectly.")
    elif mode == 3:
        if len(line) <= 6:
            warning(".ped missing genotypes.")
            
        # Check the map file as well
        with open(filename.split('.')[0] + ".map", 'r') as infile:
            line = infile.readline().split()
        if len(line) != 4:
            warning(".map incorrectly formatted.")
    else:
        warning("Unrecognized file format.")
    
    print("[ ", filename, " ] is appropriately formatted.")

###############################################################################
def getStats(filename):
    print("Reading [ ", filename, " ].")
    
    stats = {}
    with open(filename, 'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.split()
            stats[line[0]] = [float(line[5]), float(line[6])]

    return stats

###############################################################################
def filterDsf(inname, outname, stats, miss, maf):
    print("Filtering [ ", inname, " ].")

    infile = open(inname, 'r')
    keepfile = open(outname + ".dsf", 'w')
    filtfile = open(outname + "_filtered.dsf", 'w')

    header = infile.readline().split()
    keepfile.write('\t'.join(header) + '\n')
    filtfile.write('\t'.join(header) + '\n')
    
    kept = filt = 0
    for snp in infile:
        snp = snp.split()
        if snp[0] not in stats:
            warning(snp[0] + " is not present in .stat file.")

        # Filter or keep
        if float(snp[3]) <= miss and float(snp[4]) >= maf:
            keepfile.write('\t'.join(snp) + '\n')
            kept += 1
        else:
            filtfile.write('\t'.join(snp) + '\n')
            filt += 1

    infile.close()
    keepfile.close()
    filtfile.close()

    print("Kept [ ", str(kept), " ] SNPs in [ ", outname + ".dsf", " ].")
    print("Removed [ ", str(filt), " ] SNPs to [ ", outname + "_filtered.dsf", " ].")

###############################################################################
if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())
    
    # Change the working directory if necessary
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input'] is None:
        warning("No input file.")
    if args['output'] is None:
        warning("No output file.")
    
    print(version())
    
    st = timeit.default_timer()
    
    # Check input file
    checkFile(args['input'], args['mode'])

    stats = getStats(args['stat'])

    if args['mode'] == 1:
        filterDsf(args['input'], args['output'], stats, args['miss'], args['maf'])
    else:
        warning("Unrecognized input mode.")
    
    et = timeit.default_timer()

    print("Conversion finished.")
    print("Time: %.2f min." % ((et - st)/60))
