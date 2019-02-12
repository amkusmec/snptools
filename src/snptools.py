"""
Global dictionaries for nucleotide codes and helper functions for snptools.
"""

import sys

iupac = { 'AT': 'W', 'TA': 'W', 'CG': 'S', 'GC': 'S',
          'AC': 'M', 'CA': 'M', 'GT': 'K', 'TG': 'K',
          'AG': 'R', 'GA': 'R', 'CT': 'Y', 'TC': 'Y',
          '+-': '0', '-+': '0', 'AA': 'A', 'CC': 'C',
          'GG': 'G', 'TT': 'T', 'NN': 'N' }

iupac2 = { 'A': 'A A', 'C': 'C C', 'G': 'G G', 'T': 'T T',
           'M': 'A C', 'R': 'A G', 'W': 'A T', 'S': 'C G',
           'Y': 'C T', 'K': 'G T', 'N': 'N N', '+': '+ +',
           '-': '- -', '0': '+ - ' }

iupac3 = { 'A T': 'W', 'T A': 'W', 'C G': 'S', 'G C': 'S',
           'A C': 'M', 'C A': 'M', 'G T': 'K', 'T G': 'K',
           'A G': 'R', 'G A': 'R', 'C T': 'Y', 'T C': 'Y',
           '+ -': '0', '- +': '0', 'A A': 'A', 'C C': 'C',
           'G G': 'G', 'T T': 'T', 'N N': 'N' }

###############################################################################
def warning(*objs):
    print("WARNING: ", *objs, end = '\n', file = sys.stderr)
    sys.exit()

###############################################################################
def checkFile(filename, modei):
    print("Checking [", filename, "].")
    with open(filename, 'r') as infile:
        line = infile.readline().split()
        
    if modei == 1:
        if line[0] != 'snpid' or line[1] != 'major' or line[2] != 'minor':
            warning("DSF formatted incorrectly.")
    elif modei == 2:
        if line[0] != 'rs' and line[0] != 'rs#':
            warning("HapMap formatted incorrectly.")
    elif modei == 3:
        if len(line) <= 6:
            warning("PED missing genotypes.")
        with open(filename[:-3] + 'map', 'r') as infile:
            line = infile.readline().split()
            if len(line) != 4:
                warning("Map file formatted incorrectly.")
    elif modei == 4:
        print("NB: No format checking for VCF files.")
    else:
        warning("Unrecognized input mode " + str(modei) + ".")
    
    print("[", filename, "] is appropriately formatted.")
