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
    (c) 2015 Aaron Kusmec
    
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
    parser.add_argument('-o', '--output', help = 'Output file stem', type = str)
    parser.add_argument('-mi', '--modei', help = 'Input mode', type  = int)
    parser.add_argument('-c', '--coding', help = '1-based if True, 0-based if False', \
                        action = "store_true")
    
    return parser

###############################################################################
def checkFile(filename, modei):
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
    
###############################################################################
def readFile(filename, modei):
    snps = []
    
    with open(filename, 'r') as infile:
        if modei == 1:
            snps.append(["taxa"] + infile.readline().split()[5:])
            
            for line in infile:
                line = line.split()
                snps.append([line[0]] + line[5:])
            
        elif modei == 2:
            snps.append(["taxa"] + infile.readline().split()[11:])
            
            for line in infile:
                line = line.split()
                snps.append([line[0]] + line[11:])
            
        elif modei == 3:
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
def numericalize(snps, coding):
    counter = 0
    num = [snps[0]]
    
    if coding:
        major = '0'
        minor = '2'
        hetero = 0
    else:
        major = '-1'
        minor = '1'
        hetero = -1

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
            maf = count2/(count1 + count2)
            x = x.replace(allele1, major)
            x = x.replace(allele2, minor)
        else:
            maf = count1/(count1 + count2)
            x = x.replace(allele1, minor)
            x = x.replace(allele2, major)
        
        x = x.replace(het, str(hetero + 1))
        x = x.replace('N', str(hetero + 2*maf))
        
        num.append([s[0]] + x.split('\t'))
    
    return num

###############################################################################
def writeFile(num, filename):
    print("Transposing SNP matrix.")
    num = zip(*num)
    
    map_written = False
    with open(filename + '.xmat', 'w') as outfile:
        for n in num:
            if map_written:
                outfile.write('\t'.join(n) + '\n')
            else:
                outfile.write('\t'.join(n) + '\n')
                with open(filename + '.map', 'w') as mapfile:
                    for n2 in n[1:]:
                        mapfile.write('\t'.join([n2.split('_')[0], n2, '0', n2.split('_')[1]]) + '\n')
                map_written = True

if __name__ == "__main__":
    parser = getParser()
    args = vars(parser.parse_args())
    
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input'] is None or args['output'] is None:
        warning("Input/output file required.")
    
    print(version())
    
    st = timeit.default_timer()
    
    checkFile(args['input'], args['modei'])
    print("Reading from [ ", args['input'], " ].")
    snps = readFile(args['input'], args['modei'])
    
    num = numericalize(snps, args['coding'])
    
    writeFile(num, args['output'])    
    
    et = timeit.default_timer()
    
    print("Analysis finished.")
    print("Time: %.2f min." % ((et - st)/60))
    
