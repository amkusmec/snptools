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
   snpstat
   (c) 2015 Aaron Kusmec
   
   Calculate missing rates and minor allele frequencies.
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
def dsfStat(filename):
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf"]]
    counter = 0
    
    with open(filename, 'r') as infile:
        # Skip the header line
        header = infile.readline()
        
        for line in infile:
            line = line.split()
            substat = [line[0], line[0].split('_')[0], line[0].split('_')[1]]
            miss = line[5:].count("N")/len(line[5:])
            
            # Get the alleles
            geno = set(line[5:]); geno = list(geno)
            for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
                if k in geno:
                    geno += iupac2[k].split()
            geno = set(geno)
            geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
            geno = list(geno)
            
            # Make sure we have two items in the allele list even if they are
            # both missing
            if len(geno) == 0:
                geno.extend(['N', 'N'])
            elif len(geno) == 1:
                geno.append('N')
            
            # Get the alleles
            allele1 = geno[0]
            allele2 = geno[1]
            
            # Assign the heterozygous genotype
            # If only one allele is present, we assign a nonsensical heterozygous
            # genotype. If only heterozygotes were present in the SNP file, the
            # code accounts for that and decomposes it into the diploid genotypes.
            if allele1 == 'N' or allele2 == 'N':
                het = 'ZZ'
            else:
                het = iupac[allele1 + allele2]
            
            # Calculate the maf
            count1 = 2*line[5:].count(allele1) + line[5:].count(het)
            count2 = 2*line[5:].count(allele2) + line[5:].count(het)
            
            if count1 >= count2 and count1 != 0:
                maf = count2/(count1 + count2)
                substat.extend([allele1, allele2, str(miss), str(maf)])
            elif count2 > count1:
                maf = count1/(count1 + count2)
                substat.extend([allele2, allele1, str(miss), str(maf)])
            else:
                # In this case everything is missing or both alleles are
                # missing and we can't calculate a maf
                maf = 0.0
                substat.extend([allele1, allele2, str(miss), str(maf)])
            
            stats.append(substat)
            
            counter += 1
            if counter % 1e5 == 0:
                print("Processed [ ", str(counter), " ] SNPs.")
    
    return stats

###############################################################################
def hmpStat(filename):
    print("Calculating statistics.")
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf"]]
    counter = 0
    
    with open(filename, 'r') as infile:
        # Skip the header line
        header = infile.readline()
        
        for line in infile:
            line = line.strip().split('\t')
            substat = [line[0], line[2], line[3]]
            miss = line[11:].count('N')/len(line[11:])
            
            # Get the alleles
            geno = set(line[11:]); geno = list(geno)
            for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
                if k in geno:
                    geno += iupac2[k].split()
            geno = set(geno)
            geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
            geno = list(geno)
            
            # Make sure we have two items in the allele list even if they are
            # both missing
            if len(geno) == 0:
                geno.extend(['N', 'N'])
            elif len(geno) == 1:
                geno.append('N')
            
            # Get the alleles
            allele1 = geno[0]
            allele2 = geno[1]
            
            # Assign the heterozygous genotype
            # If only one allele is present, we assign a nonsensical heterozygous
            # genotype. If only heterozygotes were present in the SNP file, the
            # code accounts for that and decomposes it into the diploid genotypes.
            if allele1 == 'N' or allele2 == 'N':
                het = 'ZZ'
            else:
                het = iupac[allele1 + allele2]
            
            # Calculate the maf
            count1 = 2*line[11:].count(allele1) + line[11:].count(het)
            count2 = 2*line[11:].count(allele2) + line[11:].count(het)
            
            if count1 >= count2 and count1 != 0:
                maf = count2/(count1 + count2)
                substat.extend([allele1, allele2, str(miss), str(maf)])
            elif count2 > count1:
                maf = count1/(count1 + count2)
                substat.extend([allele2, allele1, str(miss), str(maf)])
            else:
                # In this case everything is missing or both alleles are
                # missing and we can't calculate a maf
                maf = 0.0
                substat.extend([allele1, allele2, str(miss), str(maf)])
            
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
    stats = [["snpid", "chr", "pos", "major", "minor", "miss", "maf"]]
    counter = 0      
    for s, m in zip(snps, smap):
        substat = [m[1], m[0], m[3]]
        miss = s.count("N")/len(s)
        
        # Get the alleles
        geno = set(s); geno = list(geno)
        for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
            if k in geno:
                geno += iupac2[k].split()
        geno = set(geno)
        geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
        geno = list(geno)
        
        # Make sure we have two items in the allele list even if they are
        # both missing
        if len(geno) == 0:
            geno.extend(['N', 'N'])
        elif len(geno) == 1:
            geno.append('N')
        
        # Get the alleles
        allele1 = geno[0]
        allele2 = geno[1]
        
        # Assign the heterozygous genotype
        # If only one allele is present, we assign a nonsensical heterozygous
        # genotype. If only heterozygotes were present in the SNP file, the
        # code accounts for that and decomposes it into the diploid genotypes.
        if allele1 == 'N' or allele2 == 'N':
            het = 'ZZ'
        else:
            het = iupac[allele1 + allele2]
        
        # Calculate the maf
        count1 = 2*s.count(allele1) + s.count(het)
        count2 = 2*s.count(allele2) + s.count(het)
        
        if count1 >= count2 and count1 != 0:
            maf = count2/(count1 + count2)
            substat.extend([allele1, allele2, str(miss), str(maf)])
        elif count2 > count1:
            maf = count1/(count1 + count2)
            substat.extend([allele2, allele1, str(miss), str(maf)])
        else:
            # In this case everything is missing or both alleles are
            # missing and we can't calculate a maf
            maf = 0.0
            substat.extend([allele1, allele2, str(miss), str(maf)])
        
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