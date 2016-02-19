# snptools

__`snptools` requires Python 3.4+__

`snptools` is a collection of command-line tools written in Python for manipulating high-density genotyping data. To test the tools, run the shell scripts in the `test` directory. Three file formats are supported: .dsf, .hmp.txt, and .ped. For information on each format, see `FORMATS`.

Explanations of all options for each tool can be obtained using

```
$ python3 tool_name.py --help
```

For tools requiring input/output modes, the following modes are recognized:

1. .dsf
2. .hmp.txt
3. .ped

Use of an undefined mode will generate an error.

### `convert.py`

```
$ python3 convert.py -p ../example -i example.hmp.txt -o test.dsf -mi 2 -mo 1
```

`convert` converts between different file formats. Be sure that you have specified the correct input and output modes. Specifying the wrong modes will generate undefined behavior.

### `numericalize.py`

```
$ python3 numericalize.py -p ../example -i example.hmp.txt -o test --additive --coding
```

`numericalize` converts a file to numeric genotypes, counting minor alleles. The tool also numerically imputes missing genotypes using two times the minor allele frequency (expected value of a binomial random variable).

### `snpstat.py`

```
$ python3 snpstat.py -p ../example -i example.dsf -o example.stat -mi 1
```

`snpstat` calculates missing rates and minor allele frequencies for a SNP file. This information is required by `filter`.

