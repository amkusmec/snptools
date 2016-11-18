# snptools

### Contents
- [Overview](#overview)
  - [File Formats](#file-formats)
  - [Getting Help](#getting-help)
- [Convert](#convert)
- [Filter](#filter)
- [Find Alleles](#find-alleles)
- [Numericalize](#numericalize)
- [SNP Statistics](#snp-statistics)

## Overview

`snptools` is a collection of command-line tools written in Python for manipulating text-formatted high-density genotyping data. Each function performs a single task that corresponds to the name of the function. `snptools` requires **Python 3.4+** but no modules beyond those that install by default.

### File Formats

`snptools` accepts three text-based file formats as input and can output in four different text-based file formats. For more information, see [`FORMATS`](blob/master/FORMATS.md). Most functions require input/output mode flags for unambiguous determination of file formats. These codes are

1. Density SNP Format (.dsf)
2. HapMap (.hmp.txt)
3. PLINK (.ped/.map)

Use of an undefined mode will generate an error.

### Getting Help

For detailed help, please refer to the sections for each function. Command-line arguments can be accessed for each function using
```
$ python3 src/function_name.py --help
```

## Convert

This function converts between file formats. Specifying the same input and output modes will generate an error. Note that `convert` will try to fill in metadata columns for the output file format *if* the input file also contains that information. Otherwise, appropriate missing values will be entered.

```
-p/--path   = [OPTIONAL] Path of the input file (default: .)
-i/--input  = Input file
-o/--ouput  = Output file
-mi/--modei = Integer format code for the input file
-mo/--modeo = Integer format code for the output file
```

## Filter

Filter removes SNPs based on several customizable criteria: minimum minor allele frequency, maximum missing rate, maximum heterozygosity, and inclusion in a list of SNPs to include. Filtering based on minor allele frequency, missing rate, and heterozygosity act as a series of AND filters: A SNP must meet all three criteria for inclusion.

Notes:
- A statistic file output by `snpstat` is required for filtering.
- Do not include a file extension on the output file stem. Output will be in the same format as the input. A file of removed SNPs will also be output with `_filtered` added to the stem.
- A list of SNPs to include should include the identifier for each SNP on a separate line.
- Retaining SNPs through a list is currently only supported for HapMap files.
- Retaining SNPs through a SNP list cannot currently be combined with other types of filters.

```
-p/--path   = [OPTIONAL] Path of the input file (default: .)
-s/--stat   = Statistic file
-i/--input  = Input file
-o/--output = Output file stem
-mi/--modei = Integer format code for the input file
-n/--miss   = [OPTIONAL] Floating point-valued maximum missing rate (default: 1.0)
-f/--maf    = [OPTIONAL] Floating point-valued minimum minor allele frequency (default: 0.0)
-ht/--het   = [OPTIONAL] Floating point-valued maximum heterozygosity (default: 1.0)
-r/--retain = [OPTIONAL] File of SNPs to retain (default: None)
```

## Find Alleles

This function identifies the major (most common) and minor (least common) alleles for each SNP. This is useful for obtaining allele information in DSF or HapMap files after converting from PLINK. If a SNP is not biallelic, the function will return N/N (both missing) for that SNP. This function is *not* defined for PLINK formatted files.

```
-p/--path   = [OPTIONAL] Path of the input file (default: .)
-i/--input  = Input file
-o/--output = Output file
-mi/--modei = Integer format code for the input file
```

## Numericalize

This function converts text-based SNP data into numerical genotype calls. For more information on the conversion process see [`FORMATS`](blob/FORMATS.md).

Notes:
- `numericalize` requires a statistic file output by `snpstat` to identify major and minor alleles.
- Do not include a file extension in the output argument. `numericalize` outputs two file which will be given appropriate extensions.

```
-p/--path      = [OPTIONAL] Path of the input file (default: .)
-s/--stat      = Statistic file
-i/--input     = Input file
-o/--output    = Output file stem
-mi/--modei    = Integer format code for the input file
-d/--dominance = [OPTIONAL] Flag to enable dominance coding (default: False)
```

## SNP Statistics

This function identifies the major and minor alleles, minor allele frequencies, missing rates, and heterozygosities for all SNPs. Completely missing, monomorphic, and SNPs with more than two alleles will be given the following missing values:
- minor allele frequency = -9
- missing rate = 9
- heterozygosity = 9

Such SNPs will be removed by the default parameters of `filter`.

```
-p/--path = [OPTIONAL] Path of the input file (default: .)
-i/--input = Input file
-o/--output = Output file
-mi/--modei = Integer format code for the input file
```