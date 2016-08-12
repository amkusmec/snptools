import numpy as np
import pandas as pd
import argparse
import sys
import textwrap
import timeit
import os

###############################################################################
def warning(*objs):
    print("WARNING:", *objs, end = '\n', file = sys.stderr)
    sys.exit()

###############################################################################
def getParser():
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter, description = textwrap.dedent("Placeholder text."))
    
    parser.add_argument('-i', '--input', help = 'Input file', type = str)
    parser.add_argument('-o', '--output', help = 'Output file stem', type = str)
    parser.add_argument('-mi', '--modei', help = 'Input mode', type = int)
    
    return parser

###############################################################################
def readFile(filename, modei):
    if modei == 1: # dsf file
        d = readDSF(filename)
    elif modei == 2: # hmp file
        d = readHMP(filename)
    else:
        warning("Unrecognized file format.")
    
    print("Loaded [", d.shape[0], "] SNPs from [", d.shape[1] - 1, "] individuals.")
    return d

###############################################################################
def readDSF(filename):
    with open(filename, 'r') as infile:
        line = infile.readline().strip().split('\t')
        if line[0] != "snpid" or line[1] != "major" or line[2] != "minor":
            warning("DSF file formatted incorrectly.")
    
    d = pd.read_table(filename, sep = '\t')
    d.drop(['major', 'minor', 'miss', 'maf'], axis = 1, inplace = True)
    return d

###############################################################################
def readHMP(filename):
    with open(filename, 'r') as infile:
        line = infile.readline().strip().split('\t')
        if line[0] != "rs" and line[0] != "rs#":
            warning("HMP file formatted incorrectly.")
    
    d = pd.read_table(filename, sep = '\t')
    d.drop(['alleles', 'chrom', 'pos', 'strand', 'assembly#', 'center', 'protLSID', 'assayLSID', 'panel', 'QCcode'], axis = 1, inplace = True)
    return d

###############################################################################
def numericalize(snps):
    print()
    counter = 0
    
    # Set up some constants for coding
    major = '0'
    minor = '2'
    hetero = '1'
    
    # Insert np.nan so Pandas recognizes missing data points
    snps.replace(to_replace = 'N', value = np.nan, inplace = True)
    
    for i in range(snps.shape[0]):
        # Progress updates
        counter += 1
        if counter % 1e5 == 0:
            print("Processed", counter, "of", snps.shape[1] - 1, "markers.", end = '\r')
        
        # Get allele counts
        allele_counts = snps.iloc[i, 1:].value_counts()
        major_allele, minor_allele, het_allele = "N", "N", "N"
        major_count, minor_count, het_count = 0, 0, 0
        
        # Identify major and minor alleles
        for k, v in allele_counts.iteritems():
            if k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
                het_allele, het_count = k, v
            elif v > major_count:
                minor_count, major_count = major_count, v
                minor_allele, major_allele = major_allele, k
            else:
                minor_allele, minor_count = k, v
        
        # Replace alleles with numerical values
        snps.iloc[i, 1:].replace(to_replace = [major_allele, minor_allele, het_allele], value = [major, minor, hetero], inplace = True)
        
    # Replace missing values with the mean of the genotypes
    for i in range(1, snps.shape[1]):
        snps.iloc[:, i] = pd.to_numeric(snps.iloc[:, i])
    avg_alleles = snps.mean(axis = 1)
    for i in range(snps.shape[0]):
        snps.iloc[i, 1:] = snps.iloc[i, 1:].fillna(avg_alleles[i])

if __name__ == "__main__":
    parser = getParser()
    args = vars(parser.parse_args())
    
    if args['input'] is None or args['output'] is None:
        warning("Input/output file required.")
    
    print("Placeholder text.")
    
    st = timeit.default_timer()
    
    print("Reading from [", args['input'], "].")
    snps = readFile(args['input'], args['modei'])
    
    print("Converting SNPs.")
    numericalize(snps)
    
    # Write the transposed SNP matrix to file
    print("\nWriting SNPs to [", args['output'], ".xmat].")
    snps.columns = ['taxa'] + list(snps.columns)[1:]
    snps.T.to_csv(args['output'] + '.xmat', sep = '\t', header = False, index = True)
    
    # Write the map file
    print("Writing map file to [", args['output'], ".map].")
    snp_ids = snps['taxa']
    with open(args['output'] + '.map', 'w') as outfile:
        outfile.write('\t'.join(['Chromosome', 'SNP', 'Genetic', 'Position']) + '\n')
        for s in snp_ids:
            outfile.write('\t'.join([s.split('_')[0], s, '0', s.split('_')[1]]) + '\n')
    
    et = timeit.default_timer()
    
    print("Job finished.")
    print("Time: %.2f min." % ((et - st)/60))
