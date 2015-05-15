#!/bin/bash

### test_convert.sh
### Test SNP file conversions.

cd ..

### .dsf -> .hmp.txt
python3 convert.py -i example/example.dsf -o example/test.hmp.txt -mi 1 -mo 2
diff example/example.hmp.txt example/test.hmp.txt

### .dsf -> .ped/.map
python3 convert.py -i example/example.dsf -o example/test.ped -mi 1 -mo 3
diff example/example.ped example/test.ped
diff example/example.map example/test.map

### `diff` command should show 4 different lines--one for each SNP--where missing
### rates and minor allele frequencies are missing due to .hmp.txt conversion.
### .hmp.txt -> .dsf
python3 convert.py -i example/example.hmp.txt -o example/test.dsf -mi 2 -mo 1
diff example/example.dsf example/test.dsf

### .hmp.txt -> .ped/.map
python3 convert.py -i example/example.hmp.txt -o example/test.ped -mi 2 -mo 3
diff example/example.ped example/test.ped
diff example/example.map example/test.map

### `diff` commands should show 4 different lines--one for each SNP--where alleles,
### missing rates, and minor allele frequencies are missing due to .ped/.map conversion.
### .ped/.map -> .dsf
python3 convert.py -i example/example.ped -o example/test.dsf -mi 3 -mo 1
diff example/example.dsf example/test.dsf

### .ped/.map -> .hmp.txt
python3 convert.py -i example/example.ped -o example/test.hmp.txt -mi 3 -mo 2
diff example/example.hmp.txt example/test.hmp.txt
