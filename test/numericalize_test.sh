#!/bin/bash

### numericalize_test.sh
### Test SNP file numericalization.

cd ..

### .dsf -> .xmat
python3 numericalize.py -i example/example.dsf -o example/test.xmat
diff example/example.xmat example/test.xmat

### .hmp.txt -> .xmat
python3 numericalize.py -i example/example.hmp.txt -o example/test.xmat
diff example/example.xmat example/test.xmat

### .ped/.map -> .xmat
python3 numericalize.py -i example/example.ped -o example/test.xmat
diff example/example.xmat example/test.xmat

