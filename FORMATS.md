# `snptools` Formats

## Contents
   - [IUPAC Nucleotide Ambiguity Codes](#iupac-nucleotide-ambiguity-codes)
   - [HapMap](#hapmap)
   - [PLINK](#plink)
   - [Density SNP Format](#density-snp-format)
   - [Numerical](#numerical)

### IUPAC Nucleotide Ambiguity Codes

`snptools` uses single-bit IUPAC nucleotide ambiguity codes to represent SNP genotypes. The standard codes are

| Single-bit Code | Diploid Genotype | Notes                             |
|:---------------:|:----------------:|:---------------------------------:|
| N               | N/N              |                                   |
| A               | A/A              |                                   |
| C               | C/C              |                                   |
| G               | G/G              |                                   |
| T               | T/T              |                                   |
| +               | +/+              | small insertion                   |
| -               | -/-              | small deletion                    |
| W               | A/T              |                                   |
| S               | C/G              |                                   |
| M               | A/C              |                                   |
| K               | G/T              |                                   |
| R               | A/G              |                                   |
| Y               | C/T              |                                   |
| 0               | +/-              | (heterozygous insertion/deletion) |

### HapMap

**File extension:** `.hmp.txt`

All SNPs are coded using single-bit ambiguity codes. `snptools` does *not* support two-bit coding. Column headers and definitions are
1. `rs` OR `rs#` - alphanumeric SNP identifier
2. `alleles` - separate by '/'; if possible `snptools` uses major/minor ordering
3. `chrom` - chromosome number
4. `pos` - physical (bp) position
5. `strand` - one of +, -, *
6. `assembly`
7. `center`
8. `protLSID` - protocol ID
9. `assayLSID`
10. `panel`
11. `QCcode`

Columns 12+ are named by sample and contain the SNP genotype calls for each sample.

### PLINK

**File extension:** `.ped` and `.map`

These formats are compatible with Sean Purcell's PLINK software <pngu.mgh.harvard.edu/~purcell/plink/index.shtml>. Both the pedigree and map files must be present in the same directory and have the same file stem in order to be recognized by `snptools`.

#### Pedigree file

All SNPs are coded using double-bit IUPAC nucleotide codes. Note that while PLINK accepts any biallelic coding scheme and arbitrary missing genotypes, `snptools` only operates on standard IUPAC codes. Columns do *not* have headers. Column definitions are
1. Family ID
2. Individual ID
3. Paternal ID
4. Maternal ID
5. Sex
6. Phenotype

Columns 7+ are biallelic SNP genotypes. Every two columns contains the SNP calls for a single site. Note that of the first six columns, `snptools` only operates on the family ID.

#### Map file

Each row contains the map information for a single SNP, presented in the same order in which they occur in the pedigree file. Columns do *not* have headers. Columns definitions are
1. Chromosome
2. Unique SNP identifier
3. Genetic distance (Morgans)
4. Physical position (bp)

### Density SNP Format

**File extension:** `.dsf`

All SNPs are coded using single-bit IUPAC codes. Columns headers and definitions are
1. `snpid` - alphanumeric identifier; typicallly chr_pos
2. `major` - major (most common) allele
3. `minor` - minor (least common) allele
4. `miss` - missing rate
5. `maf` - minor allele frequency

Columns 6+ are named by sample and contain the SNP calls for each sample.

### Numerical

**File extension:** `.xmat` and `.map`

This format requires both the matrix and map files to be in the same directory and have the same file stem. The map file uses the same format as the PLINK map file.

The matrix file has an initial column with header `taxa` that contains the sample IDs. Every other column has a header corresponding to a SNP identifier and contains the numerical genotype calls for that SNP. Matrix files can be created using additive or dominance coding.

#### Additive coding

This coding counts the number of minor alleles that a sample has. Assume alleles A and B where B is the minor allele. Permissible values are in the interval [0,2] according to
- AA = 0
- AB = 1
- BB = 2

Missing genotypes are imputed with the mean non-missing genotype for that SNP, or 2*p where p is the minor allele frequency.

#### Dominance coding

This coding is 1 if the sample is heterozygous at a SNP and 0 otherwise. Missing genotypes are imputed with the mean non-missing genotype, or p where p is the frequency of heterozygotes.
