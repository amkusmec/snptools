# -*- coding: utf-8 -*-
"""
Created on Fri May 22 08:29:35 2015

@author: aaron
"""

import sys
import argparse
import textwrap
import timeit
import os

# Define a dictionary converting IUPAC multi-base codes
iupac = { 'AT': 'W', 'TA': 'W',
          'CG': 'S', 'GC': 'S',
          'AC': 'M', 'CA': 'M',
          'GT': 'K', 'TG': 'K',
          'AG': 'R', 'GA': 'R',
          'CT': 'Y', 'TC': 'Y',
          '+-': '0', '-+': '0',
          'AA': 'A', 'CC': 'C',
          'GG': 'G', 'TT': 'T',
          'NN': 'N'
          }

iupac2 = { 'A': 'A A', 'C': 'C C', 'G': 'G G', 'T': 'T T',
           'M': 'A C', 'R': 'A G', 'W': 'A T', 'S': 'C G',
           'Y': 'C T', 'K': 'G T', 'N': 'N N', '+': '+ +',
           '-': '- -', '0': '+ -'
         }

###############################################################################
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

###############################################################################
def version():
   v0 = """
   ############################################################################
   snpstat
   (c) 2015 Aaron Kusmec
   
   Calculate missing rate, minor allele frequencie, and HWE.
   Input modes,
       1 = .dsf
       2 = .hmp.txt
       3 = .ped (PLINK)
       
   Usage: python3 snpstat.py -i example.dsf -o example.stat -mi 1
   
   ----------------------------------------------------------------------------
   
   Todo:
       - hmp and ped support
       - HWE calculations
       
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
    parser.add_argument('-i', '--input', help = 'Input file', type = str)
    parser.add_argument('-o', '--output', help = 'Output file', type = str)
    parser.add_argument('-mi', '--modei', help = 'Input mode', type = int)

    return parser

###############################################################################
def checkFile(filename, modei):
    print("Checking [ ", filename, " ].")
    
    with open(filename, 'r') as infile:
        line = infile.readline().split()
    
    if modei == 1:
        if line[0] != "snpid" or line[1] != "major" or line[2] != "minor":
            warning(".dsf formatted incorrectly.")
    elif modei == 2:
        if line[0] != "rs" and line[0] != "rs#":
            warning(".hmp.txt formatted incorrectly.")
    elif modei == 3:
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
def writeStats(stats, filename):
    print("Writing [ ", filename, " ].")
    with open(filename, 'w') as outfile:
        for s in stats:
            outfile.write('\t'.join(s) + '\n')

###############################################################################
def dsfStat(filename):
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf"]]
    counter = 0
    
    with open(filename, 'r') as infile:
        # Skip the header line
        header = infile.readline()
        
        for line in infile:
            counter += 1
            if counter % 1e5 == 0:
                print("Processed [ ", str(counter), " ] SNPs.")
            
            line = line.split()
            substat = [line[0], line[0].split('_')[0], line[0].split('_')[1]]
            miss = line[5:].count("N")/len(line[5:])
            
            # Get the alleles and maf
            geno = set(line[5:])
            geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
            geno = list(geno)
            allele1 = geno[0]
            allele2 = geno[1]
            het = iupac[allele1 + allele2]
            
            count1 = 2*line[5:].count(allele1) + line[5:].count(het)
            count2 = 2*line[5:].count(allele2) + line[5:].count(het)
            
            if count1 >= count2:
                maf = count2/(count1 + count2)
                substat.extend([allele1, allele2, str(miss), str(maf)])
            else:
                maf = count1/(count1 + count2)
                substat.extend([allele2, allele1, str(miss), str(maf)])
            
            stats.append(substat)
    
    return stats

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
    checkFile(args['input'], args['modei'])
    
    if args['modei'] == 1:
        stats = dsfStat(args['input'])
    else:
        warning("Unrecognized input mode.")
        
    # Write output
    writeStats(stats, args['output'])
    
    et = timeit.default_timer()

    print("Conversion finished.")
    print("Time: %.2f min." % ((et - st)/60))