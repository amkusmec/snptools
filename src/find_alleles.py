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
def findAllelesDsf(inname, outname):
    with open(inname, 'r') as infile:
        with open(outname, 'w') as outfile:
            header = infile.readline()
            outfile.write(header)

            for line in infile:
                line = line.split('\t')
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
                line = line.split('\t')
                major, minor = majorMinor(line[11:])
                line[1] = major + '/' + minor
                outfile.write('\t'.join(line) + '\n')

###############################################################################
def majorMinor(snp):
    set0 = set(snps)
    alleles = list(set0)

    for k in ['W', 'S', 'M', 'K', 'R', 'Y', '0']:
        if k in alleles:
            alleles += iupac2[k].split()

    set0 = set(alleles)
    set0 = set0.difference(set(['W', 'S', 'M', 'K', 'R', 'Y', '0', 'N']))
    alleles = list(set0)

    if len(alleles) == 2:
        het = iupac[''.join(alleles)]
        allele1 = 2*snps.count(alleles[0]) + snps.count(het)
        allele2 = 2*snps.count(alleles[1]) + snps.count(het)

        if allele1 >= allele2:
            return alleles[0], alleles[1]
        else:
            return alleles[1], alleles[0]

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
    snps = readFile(args['input'], args['modei'])
    
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








