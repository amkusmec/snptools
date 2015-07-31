import sys
import argparse
import textwrap
import timeit
import os
from operator import itemgetter

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
   sort
   (c) 2015 Aaron Kusmec
   
   Sorts the input file by chromosome and base-pair position.
   
   Limitations:
    - Only works for .dsf files
    - Output only in the input format
   
   Usage: python3 sort.py
   
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
def sortDsf(snps):
    # Remove the header line
    header = snps[0]
    snps = snps[1:]

    # Create the key and place at the end of each line
    for i in range(len(snps)):
        key = [int(snps[i][0].split('_')[0]), int(snps[i][0].split('_')[1])]
        snps[i].extend(key)

    # Use Python's built-in list sort to sort the SNPs in place
    snps.sort(key = itemgetter(-2, -1))

    # Remove the key
    for i in range(len(snps)):
        snps[i] = snps[i][:-2]

    snps.insert(0, header)
    return snps

###############################################################################
def sortHmp(snps):
    # Remove the header line
    header = snps[0]
    snps = snps[1:]

    # Split chromosome and position numbers already exist
    # We just need to convert them to integers
    for i in range(len(snps)):
        snps[i][2] = int(snps[i][2])
        snps[i][3] = int(snps[i][3])

    # Use Python's built-in list sort to sort the SNPs in place
    snps.sort(key = itemgetter(2, 3))

    # Convert back to strings
    for i in range(len(snps)):
        snps[i][2] = str(snps[i][2])
        snps[i][3] = str(snps[i][3])

    snps.insert(0, header)
    return snps

###############################################################################
if __name__ == '__main__':
    parser = get_parser()
    args = vars(parser.parse_args())
    
    # Change the working directory if necessary
    if args['path'] is not None:
        os.chdir(args['path'])
    if args['input'] is None:
        warning("Missing input file(s).")
    if args['output'] is None:
        warning("No output file.")
    
    print(version())
    
    st = timeit.default_timer()
    
    # Load files
    checkFile(args['input'], args['modei'])
    snps = readFile(args['input'], args['modei'])
    
    # Merge
    if args['modei'] == 1:
        sorted = sortDsf(snps)
    if args['modei'] == 2:
        sorted = sortHmp(snps)
    else:
        warning("Unrecognized input/output mode(s).")
        
    # Write output
    writeFile(sorted, args['output'])
    
    et = timeit.default_timer()

    print("Sort finished.")
    print("Time: %.2f min." % ((et - st)/60))
    















