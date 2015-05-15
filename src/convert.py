# -*- coding: utf-8 -*-
"""
Created on Fri May 15 11:17:44 2015

@author: Aaron Kusmec
"""

import sys
import argparse
import textwrap
import timeit
import os

# Define a dictionary converting IUPAC multi-base codes
iupac = { 'A T': 'W', 'T A': 'W',
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
   convert
   (c) 2015 by Aaron Kusmec
   
   Convert between different SNP file formats.
   Input/output modes,
       1 = .dsf
       2 = .hmp.txt
       3 = .ped (PLINK)
       
   Usage: python3 convert.py -i example.dsf -o example.ped -mi 1 -mo 3
   
   ----------------------------------------------------------------------------
   
   NOTE
   When converting between file types, metadata fields that do not exist in 
   both file types will be replaced with appropriate missing values. If you 
   need allele information in .dsf or .hmp.txt files, run the `find_alleles` 
   utility on the output. If you need missing rates and minor allele 
   frequencies in .dsf files, run the `` utility on the output.
   
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
    parser.add_argument('-mo', '--modeo', help = 'Output mode', type = int)

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
def readFile(filename, modei):
    snps = []
    
    print("Reading from [ ", filename, " ].")
    with open(filename, 'r') as infile:
        for line in infile:
            snps.append(line.strip().split('\t'))
    
    if modei == 3:
        smap = []
        
        print("Reading from [ ", filename.split('.')[0] + ".map", " ].")
        with open(filename.split('.')[0] + ".map", 'r') as infile:
            for line in infile:
                smap.append(line.split())
    else:
        smap = None
    
    return snps, smap

###############################################################################
def writeFile(snps, smap, filename):
    if smap is not None:
        print("Writing [ ", filename.split('.')[0] + ".map", " ].")
        with open(filename.split('.')[0] + ".map", 'w') as outfile:
            for s in smap:
                outfile.write('\t'.join(s) + '\n')
    
    print("Writing [ ", filename, " ].")
    with open(filename, 'w') as outfile:
        for s in snps:
            outfile.write('\t'.join(s) + '\n')

###############################################################################
def dsf2hmp(snps):
    print("Read [ ", str(len(snps) - 1), " ] SNPs from [ ", 
                 str(len(snps[0][5:])), " ] samples.")
    print("Converting from .dsf to .hmp.txt.")
    
    snps = zip(*snps)
    
    ### Create new metadata
    # Process IDs
    snpid = next(snps)
    rs = ["rs#"] + list(snpid[1:])
    snpid = [s.split('_') for s in snpid[1:]]
    chrom = ["chrom"] + [s[0] for s in snpid]
    pos = ["pos"] + [s[1] for s in snpid]
    
    # Alleles
    major = next(snps); minor = next(snps)
    alleles = ["alleles"] + ['/'.join(z) for z in zip(major[1:], minor[1:])]
    
    # Other
    strand = ["strand"] + ["NA"]*(len(rs) - 1)
    assembly = ["assembly#"] + ["NA"]*(len(rs) - 1)
    center = ["center"] + ["NA"]*(len(rs) - 1)
    prot = ["protLSID"] + ["NA"]*(len(rs) - 1)
    assay = ["assayLSID"] + ["NA"]*(len(rs) - 1)
    panel = ["panel"] + ["NA"]*(len(rs) - 1)
    qc = ["QCcode"] + ["NA"]*(len(rs) - 1)
    
    ### Remove the last metadata columns from the .dsf input
    miss = next(snps); maf = next(snps)
    
    ### Reformat the snp table
    snps = zip(rs, alleles, chrom, pos, strand, assembly, center, prot, assay,
               panel, qc, *snps)
    return snps

###############################################################################
def hmp2dsf(snps):
    print("Read [ ", str(len(snps) - 1), " ] SNPs from [ ",
                 str(len(snps[0][11:])), " ] samples.")
    print("Converting from .hmp.txt to .dsf.")
    
    snps = zip(*snps)
    
    ### Create new metadata
    rs = next(snps)
    
    # Process alleles
    alleles = next(snps)
    alleles = [a.split('/') for a in alleles[1:]]
    major = ["major"] + [a[0] for a in alleles]
    minor = ["minor"] + [a[1] for a in alleles]
    
    # Process IDs
    chrom = next(snps); pos = next(snps)
    snpid = ["snpid"] + ['_'.join(z) for z in zip(chrom[1:], pos[1:])]
    
    # Other
    miss = ["miss"] + ["0"]*len(snpid)
    maf = ["maf"] + ["0"]*len(snpid)
    
    ### Remove the last metadata
    for _ in range(4, 11):
        temp = next(snps)
    
    ### Reformat the snp table
    snps = zip(snpid, major, minor, miss, maf, *snps)
    return snps

###############################################################################
def dsf2ped(snps):
    print("Read [ ", str(len(snps) - 1), " ] SNPs from [ ",
                 str(len(snps[0][5:])), " ] samples.")
    print("Converting from .dsf to .ped/.map.")
    
    ### Get the IDs before transposition
    ids = snps[0][5:]
    
    ### Genotypes
    for i in range(1, len(snps)):
        snps[i] = snps[i][:5] + [iupac2[s] for s in snps[i][5:]]
    
    snps = zip(*snps[1:])
    
    ### Create the .map file
    snpid = next(snps)
    snpid = [s.split('_') for s in snpid]
    chrom = [s[0] for s in snpid]
    pos = [s[1] for s in snpid]
    snpid = ['_'.join(s) for s in snpid]
    gen = ["0"]*len(snpid)
    
    smap = zip(chrom, snpid, gen, pos)
    
    ### Remove metadata
    for _ in range(1, 5):
        temp = next(snps)
    
    ### Create the .ped file
    ## Metadata
    pid = ["0"]*len(ids)
    mid = ["0"]*len(ids)
    sex = ["0"]*len(ids)
    pheno = ["-9"]*len(ids)
    
    snps = zip(*snps)
    snps = zip(ids, ids, pid, mid, sex, pheno, *snps)
    return snps, smap

###############################################################################
def ped2dsf(snps, smap):
    print("Read [ ", str(len(smap)), " ] SNPs from [ ",
                 str(len(snps)), " ] samples.")
    print("Converting from .ped/.map to .dsf.")
    
    ### Get the snp IDs
    smap = zip(*smap)
    chrom = next(smap); rs = next(smap); gen = next(smap); pos = next(smap)
    snpid = ["snpid"] + ['_'.join(z) for z in zip(chrom, pos)]
    
    ### Genotypes
    for i in range(len(snps)):
        snps[i] = [snps[i][1]] + [iupac[s] for s in snps[i][6:]]
    
    ### Other metadata
    major = ["major"] + ["N"]*len(snps)
    minor = ["minor"] + ["N"]*len(snps)
    miss = ["miss"] + ["0"]*len(snps)
    maf = ["maf"] + ["0"]*len(snps)
    
    snps = zip(snpid, major, minor, miss, maf, *snps)
    return snps

###############################################################################
def hmp2ped(snps):
    print("Read [ ", str(len(snps) - 1), " ] SNPs from [ ",
                 str(len(snps[0][11:])), " ] samples.")
    print("Converting from .hmp.txt to .ped/.map.")
    
    ### Get the IDs before transposition
    ids = snps[0][11:]
    
    ### Genotypes
    for i in range(1, len(snps)):
        snps[i] = snps[i][:11] + [iupac2[s] for s in snps[i][11:]]
    
    snps = zip(*snps[1:])
    
    ### Create the .map file
    rs = next(snps); alleles = next(snps)
    chrom = next(snps); pos = next(snps)
    
    gen = ["0"]*len(rs)
    
    smap = zip(chrom, rs, gen, pos)
    
    ### Remove metadata
    for _ in range(4, 11):
        temp = next(snps)
    
    ### Create the .ped file
    ## Metadata
    pid = ["0"]*len(ids)
    mid = ["0"]*len(ids)
    sex = ["0"]*len(ids)
    pheno = ["-9"]*len(ids)
    
    snps = zip(*snps)
    snps = zip(ids, ids, pid, mid, sex, pheno, *snps)
    return snps, smap

###############################################################################
def ped2hmp(snps, smap):
    print("Read [ ", str(len(smap)), " ] SNPs from [ ",
                 str(len(snps)), " ] samples.")
    print("Converting from .ped/.map to .hmp.txt.")
    
    ### Create new metadata
    smap = zip(*smap)
    chrom = ["chrom"] + list(next(smap))
    rs = ["rs#"] + list(next(smap))
    gen = next(smap)
    pos = ["pos"] + list(next(smap))
    alleles = ["alleles"] + ["N/N"]*(len(rs) - 1)
    strand = ["strand"] + ["NA"]*(len(rs) - 1)
    assembly = ["assembly#"] + ["NA"]*(len(rs) - 1)
    center = ["center"] + ["NA"]*(len(rs) - 1)
    prot = ["protLSID"] + ["NA"]*(len(rs) - 1)
    assay = ["assayLSID"] + ["NA"]*(len(rs) - 1)
    panel = ["panel"] + ["NA"]*(len(rs) - 1)
    qc = ["QCcode"] + ["NA"]*(len(rs) - 1)
    
    ### Genotypes
    for i in range(len(snps)):
        snps[i] = [snps[i][1]] + [iupac[s] for s in snps[i][6:]]
    
    snps = zip(rs, alleles, chrom, pos, strand, assembly, center,
               prot, assay, panel, qc, *snps)
    return snps

###############################################################################
if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())
    
    # Change the working directory if necessary
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input'] is None:
        warning("No input file!")
    if args['output'] is None:
        warning("No output file!")

    st = timeit.default_timer()
    
    # Load file
    checkFile(args['input'], args['modei'])
    snps, smap = readFile(args['input'], args['modei'])
    
    # Convert
    if args['modei'] == 1 and args['modeo'] == 2:
        snps = dsf2hmp(snps)
    elif args['modei'] == 2 and args['modeo'] == 1:
        snps = hmp2dsf(snps)
    elif args['modei'] == 1 and args['modeo'] == 3:
        snps, smap = dsf2ped(snps)
    elif args['modei'] == 3 and args['modeo'] == 1:
        snps = ped2dsf(snps, smap)
        smap = None
    elif args['modei'] == 2 and args['modeo'] == 3:
        snps, smap = hmp2ped(snps)
    elif args['modei'] == 3 and args['modeo'] == 2:
        snps = ped2hmp(snps, smap)
        smap = None
    else:
        warning("Unrecognized input/output mode(s).")
        
    # Write output
    writeFile(snps, smap, args['output'])
    
    et = timeit.default_timer()

    print("Conversion finished.")
    print("Time: %.2f min." % ((et - st)/60))