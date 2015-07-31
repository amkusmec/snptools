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
   merge
   (c) 2015 Aaron Kusmec
   
   implements a union join on taxa and snps - includes every unique taxa-snp
   combination; every combination not in a file is filled in by a missing value
   
   Limitations:
    - Only works for .dsf files
    - Both input files must be the same format
    - Output only in the input format
   
   Usage: python3 merge.py
   
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
    parser.add_argument('-i1', '--input1', help = 'Input file 1', type = str)
    parser.add_argument('-i2', '--input2', help = 'Input file 2', type = str)
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
    
    return snps

###############################################################################
def writeFile(snps, filename):
    print("Writing [ ", filename, " ].")
    with open(filename, 'w') as outfile:
        for s in snps:
            outfile.write('\t'.join(s) + '\n')

###############################################################################
def mergeDsf(snps1, snps2):
    pass

###############################################################################
if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())
    
    # Change the working directory if necessary
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input1'] is None or args['input2'] is None:
        warning("Missing input file(s).")
    if args['output'] is None:
        warning("No output file.")
    
    print(version())
    
    st = timeit.default_timer()
    
    # Load files
    checkFile(args['input1'], args['modei'])
    checkFile(args['input1'], args['modei'])
    snps1 = readFile(args['input1'], args['modei'])
    snps2 = readFile(args['input2'], args['modei'])
    
    # Merge
    if args['modei'] == 1:
        merged = mergeDsf(snps1, snps2)
    else:
        warning("Unrecognized input/output mode(s).")
        
    # Write output
    writeFile(merged, args['output'])
    
    et = timeit.default_timer()

    print("Merge finished.")
    print("Time: %.2f min." % ((et - st)/60))
    















