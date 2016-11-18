import argparse
import textwrap
import timeit
import os
from snptools import *

###############################################################################
def version():
   v0 = """
   ############################################################################
   find_alleles
   (c) 2015 Aaron Kusmec
   
   Identify major and minor alleles. Only works for .dsf and .hmp.txt files.
   Caution: Be sure to specify separate input and output files, otherwise you
        will overwrite the input file.
       
   Usage: python3 find_alleles.py -i example.dsf -o example2.dsf -mi 1
       
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
def findAllelesDsf(inname, outname):
    with open(inname, 'r') as infile:
        with open(outname, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)

            for line in infile:
                line = line.strip().split('\t')
                major, minor = majorMinor(line[5:])
                line[1], line[2] = major, minor
                outfile.write('\t'.join(line) + '\n')

###############################################################################
def findAllelesHmp(inname, outname):
    with open(inname, 'r') as infile:
        with open(outname, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)

            for line in infile:
                line = line.strip().split('\t')
                major, minor = majorMinor(line[11:])
                line[1] = major + '/' + minor
                outfile.write('\t'.join(line) + '\n')

###############################################################################
def majorMinor(snp):
    set0 = set(snp)
    alleles = list(set0)

    for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
        if k in alleles:
            alleles += iupac2[k].split()

    set0 = set(alleles)
    set0 = set0.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
    alleles = list(set0)

    if len(alleles) == 2:
        het = iupac[''.join(alleles)]
        allele1 = 2*snp.count(alleles[0]) + snp.count(het)
        allele2 = 2*snp.count(alleles[1]) + snp.count(het)

        if allele1 >= allele2:
            return alleles[0], alleles[1]
        else:
            return alleles[1], alleles[0]
    else:
        return 'N', 'N'

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
    
    # Load file
    checkFile(args['input'], args['modei'])
    
    # Convert
    if args['modei'] == 1:
        findAllelesDsf(args['input'], args['output'])
    elif args['modei'] == 2:
        findAllelesHmp(args['input'], args['output'])
    else:
        warning("Unrecognized input/output mode(s).")
    
    et = timeit.default_timer()

    print("Conversion finished.")
    print("Time: %.2f min." % ((et - st)/60))
