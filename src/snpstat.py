# -*- coding: utf-8 -*-
"""
Created on Fri May 22 08:29:35 2015

@author: aaron
"""

import argparse
import textwrap
import timeit
import os
from snptools import *

###############################################################################
def version():
   v0 = """
   ############################################################################
   snpstat - V1.5
   (c) 2015 Aaron Kusmec
   
   Calculate missing rates, minor allele frequencies, and heterozygosity.
   Completely missing, monomorphic, and triallelic+ SNPs will be given values
       of maf = -9, miss = 9, het = 9. These SNPs can be removed with filter
       and appropriate values.
   Input modes,
       1 = .dsf
       2 = .hmp.txt
       3 = .ped (PLINK)
       
   Usage: python3 snpstat.py -i example.dsf -o example.stat -mi 1
       
   ############################################################################
   """
   
   return v0

#############################################################################  
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
def writeStats(stats, filename):
    print("Writing [ ", filename, " ].")
    with open(filename, 'w') as outfile:
        for s in stats:
            outfile.write('\t'.join(s) + '\n')

###############################################################################
def calculate(snp):
    # Missing rate is simple to get
    miss = snp.count('N')/len(snp)
    
    # Find the alleles
    geno = list(set(snp))
    for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
        if k in geno:
            geno += iupac2[k].split()
    geno = set(geno)
    geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
    geno = list(geno)
    
    # Case where every site is missing
    if len(geno) == 0:
        return ['N', 'N', str(9), str(-9), str(9)]
    # Case where the SNP is monomorphic
    elif len(geno) == 1:
        return [geno[0], 'N', str(9), str(-9), str(9)]
    # Case where the SNP is triallelic or better
    elif len(geno) > 2:
        return [geno[0], ','.join(geno[1:]), str(9), str(-9), str(9)]
    
    # If we made it here, the SNP is biallelic
    allele1, allele2 = geno[0], geno[1]
    
    # Check for an insertion/deletion without proper coding
    if (allele1 == '+' or allele1 == '-') and \
       (allele2 != '+' and allele2 != '-'):
       if allele1 == '+':
          allele2 = '-'
       else:
          allele2 = '+'
    elif (allele2 == '+' or allele2 == '-') and \
         (allele1 != '+' and allele1 != '-'):
       if allele2 == '+':
          allele1 = '-'
       else:
          allele1 = '+'
    
    het = iupac[allele1 + allele2]
    count1 = 2*snp.count(allele1) + snp.count(het)
    count2 = 2*snp.count(allele2) + snp.count(het)
    hetero = snp.count(het)/(len(snp) - snp.count('N'))
    
    if count1 >= count2:
        maf = count2/(count1 + count2)
        return [allele1, allele2, str(miss), str(maf), str(hetero)]
    else:
        maf = count1/(count1 + count2)
        return [allele2, allele1, str(miss), str(maf), str(hetero)]

###############################################################################
def dsfStat(filename):
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf", "het"]]
    counter = 0
    
    with open(filename, 'r') as infile:
        # Skip the header line
        header = infile.readline()
        
        for line in infile:
            line = line.split()
            substat = [line[0], line[0].split('_')[0], line[0].split('_')[1]]
            substat.extend(calculate(line[5:]))
            stats.append(substat)
            
            counter += 1
            if counter % 1e5 == 0:
                print("Processed [ ", str(counter), " ] SNPs.")
    
    return stats

###############################################################################
def hmpStat(filename):
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf", "het"]]
    counter = 0
    
    with open(filename, 'r') as infile:
        # Skip the header line
        header = infile.readline()
        
        for line in infile:
            line = line.strip().split('\t')
            substat = [line[0], line[2], line[3]]
            substat.extend(calculate(line[11:]))
            stats.append(substat)
            
            counter += 1
            if counter % 1e5 == 0:
                print("Processed [ ", str(counter), " ] SNPs.")
    
    return stats

###############################################################################
def pedStat(filename):
    print("Reading [ ", filename.split('.')[0] + ".map", " ].")
    smap = []
    with open(filename.split('.')[0] + ".map", 'r') as infile:
        for line in infile:
            smap.append(line.split())
    
    print("Reading [ ", filename, " ].")
    snps = []
    with open(filename, 'r') as infile:
        for line in infile:
            tsnp = []
            line = line.strip().split()
            for i in range(6, len(line), 2):
               tsnp.append(iupac[line[i] + line[i+1]])
               snps.append(tsnp)
    snps = zip(*snps)
      
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf", "het"]]
    counter = 0      
    for s, m in zip(snps, smap):
        substat = [m[1], m[0], m[3]]
        substat.extend(calculate(s))
        stats.append(substat)
        
        counter += 1
        if counter % 1e5 == 0:
            print("Processed [ ", str(counter), " ] SNPs.")
    
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
    elif args['modei'] == 2:
        stats = hmpStat(args['input'])
    elif args['modei'] == 3:
        stats = pedStat(args['input'])
    else:
        warning("Unrecognized input mode.")
        
    # Write output
    writeStats(stats, args['output'])
    
    et = timeit.default_timer()

    print("Conversion finished.")
    print("Time: %.2f min." % ((et - st)/60))
