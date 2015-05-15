#!/bin/bash

### numericalize_test.sh
### Test SNP file numericalization.

cd ../src

### .dsf -> .xmat
python3 numericalize.py -p ../example -i example.dsf -o test.xmat
diff ../example/example.xmat ../example/test.xmat

### .hmp.txt -> .xmat
python3 numericalize.py -p ../example -i example.hmp.txt -o test.xmat
diff ../example/example.xmat ../example/test.xmat

### .ped/.map -> .xmat
python3 numericalize.py -p ../example -i example.ped -o test.xmat
diff ../example/example.xmat ../example/test.xmat

