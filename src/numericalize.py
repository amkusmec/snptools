# -*- coding: utf-8 -*-
"""
Created on Wed May 13 15:11:53 2015

@author: Aaron Kusmec
"""

import argparse
import textwrap
import timeit
import os
from snptools import *

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
def getParser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(version())
        )
    
    parser.add_argument('-p', '--path', help = 'Path of the input files', \
                        nargs='?', default=os.getcwd())
    parser.add_argument('-s', '--stat', help = 'SNP statistic file', type = str)
    parser.add_argument('-i', '--input', help = 'Input file', type = str)
    parser.add_argument('-o', '--output', help = 'Output file stem', type = str)
    parser.add_argument('-mi', '--modei', help = 'Input mode', type  = int)
    parser.add_argument('-d', '--dominance', help = 'Enable dominance coding', \
                        action = "store_true", default = False)
    
    return parser

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
def readStat(filename):
    print("Reading [ ", filename, " ].")

    stats = {}
    with open(filename, 'r') as infile:
        header = infile.readline()
        for line in infile:
            line = line.split()
            stats[line[0]] = { 'id' : line[0], 'chr' : line[1], 'pos' : line[2], 'major' : line[3], 'minor' : line[4], 'maf' : float(line[6]), 'het' : float(line[7]) }
    
    return stats

###############################################################################
def numericalizeA(snps, stats):
    counter = 0
    num = [snps[0]]
    
    for s in snps[1:]:
        counter += 1
        if counter % 1e5 == 0:
            print("Processed ", str(counter), " of ", str(len(snps) - 1), " markers.")
        
        if s[0] not in stats:
            warning(s[0] + " is not present in .stat file.")
        
        # Do the numerical conversion
        x = '\t'.join(s[1:])
        x = x.replace(stats[s[0]]['major'], '0')
        x = x.replace(stats[s[0]]['minor'], '2')
        
        if stats[s[0]]['major'] != 'N' and stats[s[0]]['minor'] != 'N':
            het = iupac[''.join([stats[s[0]]['major'], stats[s[0]]['minor']])]
            x = x.replace(het, '1')
        
        x = x.replace('N', str(2*stats[s[0]]['maf']))
        
        num.append([s[0]] + x.split('\t'))
    
    return num

###############################################################################
def numericalizeD(snps, stats):
    counter = 0
    num = [snps[0]]
    
    for s in snps[1:]:
        counter += 1
        if counter % 1e5 == 0:
            print("Processed ", str(counter), " of ", str(len(snps) - 1), " markers.")
        
        if s[0] not in stats:
            warning(s[0] + " is not in the .stat file.")
        
        # Get the heterozygous code
        het = iupac[''.join([stats[s[0]]['major'], stats[s[0]]['minor']])]
        
        # Perform the numerical conversion
        x = '\t'.join(s[1:])
        x = x.replace(stats[s[0]]['major'], '0')
        x = x.replace(stats[s[0]]['minor'], '0')
        x = x.replace(het, '1')
        x = x.replace('N', str(stats[s[0]]['het']))
        
        num.append([s[0]] + x.split('\t'))
        
    return num

###############################################################################
def writeFile(num, stats, filename):
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
                        mapfile.write('\t'.join([stats[n2]['chr'], stats[n2]['id'], '0', stats[n2]['pos']]) + '\n')
                    map_written = True

if __name__ == "__main__":
    parser = getParser()
    args = vars(parser.parse_args())
    
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input'] is None or args['output'] is None:
        warning("Input/output file required.")
    if args['stat'] is None:
        warning("Statistic file required.")
    
    print(version())
    
    st = timeit.default_timer()
    
    checkFile(args['input'], args['modei'])
    print("Reading from [ ", args['input'], " ].")
    snps = readFile(args['input'], args['modei'])
    
    stats = readStat(args['stat'])
    
    if args['dominance']:
        print("Using dominance coding.")
        num = numericalizeD(snps, stats)
    else:
        print("Using additive coding.")
        num = numericalizeA(snps, stats)
    
    writeFile(num, stats, args['output'])    
    
    et = timeit.default_timer()
    
    print("Analysis finished.")
    print("Time: %.2f min." % ((et - st)/60))
    
