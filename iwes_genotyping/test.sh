run -it -v /Users:/Users -v /Volumes:/Volumes dockerreg.chtc.wisc.edu/dhoconno/iwes:27329 bin/bash

cd /Volumes/T7/MHC_genotyper/baylor31/input_2
cp ~/github/iwes_genotyping/static_files/* ./
/Users/dabaker3/github/iwes_genotyping/modules/launch_iwes_smk.sh -c 2 -s 102923.fastq.gz -z T