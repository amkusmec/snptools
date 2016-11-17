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

iupac3 = { 'A T': 'W', 'T A': 'W',
           'C G': 'S', 'G C': 'S',
           'A C': 'M', 'C A': 'M',
           'G T': 'K', 'T G': 'K',
           'A G': 'R', 'G A': 'R',
           'C T': 'Y', 'T C': 'Y',
           '+ -': '0', '- +': '0',
           'A A': 'A', 'C C': 'C',
           'G G': 'G', 'T T': 'T',
           'N N': 'N'
           }

###############################################################################
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()

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
            line = line.strip().split('\t')[6:]
            snps.append([iupac3[l] for l in line])
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
