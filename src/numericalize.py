# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:11:53 2015

@author: Aaron Kusmec
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

###############################################################################
def version():
    v0 = """
    ###########################################################################
    numericalize
    (c) 2015 by Aaron Kusmec
    
    ---------------------------------------------------------------------------
    
    Numericalize a SNP genotype file in .dsf, .hmp.txt, or .ped formats. 
    Missing genotypes are replaced by 2*p where p is the frequency of the 
    minor allele (i.e., the expected value of a binomial random variable).
    
    ###########################################################################
    """
    
    return v0

###############################################################################
def warning(*objs):
    print("WARNING: ", *objs, end='\n', file=sys.stderr)
    sys.exit()
    
###############################################################################
def getParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )
    
    parser.add_argument('-p', '--path', help = 'Path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-i', '--input', help = 'Input file', type = str)
    parser.add_argument('-o', '--output', help = 'Output file', type = str)
    
    return parser

###############################################################################
def checkFile(filename):
    with open(filename, 'r') as infile:
        line = infile.readline().split()
    
    if "dsf" in filename.split('.')[1]:
        if line[0] != "snpid" or line[1] != "major" or line[2] != "minor":
            warning(".dsf formatted incorrectly.")
    elif "hmp" in filename.split('.')[1]:
        if line[0] != "rs" and line[0] != "rs#":
            warning(".hmp.txt formatted incorrectly.")
    elif "ped" in filename.split('.')[1]:
        if len(line) <= 6:
            warning(".ped missing genotypes.")
            
        # Check the map file as well
        with open(filename.split('.')[0] + ".map", 'r') as infile:
            line = infile.readline().split()
        if len(line) != 4:
            warning(".map incorrectly formatted.")
    else:
        warning("Unrecognized file format.")
    
###############################################################################
def readFile(filename):
    snps = []
    ext = filename.split('.')[1]
    
    with open(filename, 'r') as infile:
        if "dsf" in ext:
            snps.append(["taxa"] + infile.readline().split()[5:])
            
            for line in infile:
                line = line.split()
                snps.append([line[0]] + line[5:])
            
        elif "hmp" in ext:
            snps.append(["taxa"] + infile.readline().split()[11:])
            
            for line in infile:
                line = line.split()
                snps.append([line[0]] + line[11:])
            
        elif "ped" in ext:
            for line in infile:
                line = line.split()
                
                tsnp = [line[1]]
                for i in range(6, len(line), 2):
                    tsnp.append(iupac[line[i] + line[i+1]])
                snps.append(tsnp)
            
            with open(filename.split('.')[0] + ".map", 'r') as infile2:
                ids = ["taxa"]
                for line in infile2:
                    ids.append(line.split()[1])
                
                snps.insert(0, ids)
            
            snps = list(map(list, zip(*snps)))
    
    print("Loaded [ ", str(len(snps) - 1), " ] SNPs from [ ", 
                   str(len(snps[0]) - 1), " ] founders.")
    return snps

###############################################################################
def numericalize(snps):
    counter = 0
    num = [snps[0]]
    
    for s in snps[1:]:
        counter += 1
        if counter % 1e5 == 0:
            print("Processed ", str(counter), " of ", str(len(snps) - 1), " markers.")
        
        geno = set(s[1:])
        geno = geno.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
        geno = list(geno)
        allele1 = geno[0]
        allele2 = geno[1]
        het = iupac[allele1 + allele2]
        
        count1 = 2*s[1:].count(allele1) + s[1:].count(het)
        count2 = 2*s[1:].count(allele2) + s[1:].count(het)
        
        x = '\t'.join(s[1:])
        
        if count1 >= count2:
            maf = count2/(len(s[1:])*2)
            x = x.replace(allele1, '0')
            x = x.replace(allele2, '2')
        else:
            maf = count1/(len(s[1:])*2)
            x = x.replace(allele1, '2')
            x = x.replace(allele2, '0')
        
        x = x.replace(het, '1')
        x = x.replace('N', str(2*maf))
        
        num.append([s[0]] + x.split('\t'))
    
    return num

###############################################################################
def writeFile(num, filename):
    print("Transposing SNP matrix.")
    num = zip(*num)
    
    with open(filename, 'w') as outfile:
        for n in num:
            outfile.write('\t'.join(n) + '\n')

if __name__ == "__main__":
    parser = getParser()
    args = vars(parser.parse_args())
    
    if args['input'] is None or args['output'] is None:
        warning("Input/output file required.")
    
    st = timeit.default_timer()
    
    checkFile(args['input'])
    print("Reading from [ ", args['input'], " ].")
    snps = readFile(args['input'])
    
    num = numericalize(snps)
    
    writeFile(num, args['output'])    
    
    et = timeit.default_timer()
    
    print("Analysis finished.")
    print("Time: %.2f min." % ((et - st)/60))
    