bait_completed=False
while getopts s:c:z: opt ; do
   case $opt in
      s) accession=$OPTARG ;;
      c) cores=$OPTARG ;;
      z) bait_completed=$OPTARG;;
      *) usage; exit 1;;
   esac
done
echo ${cores}
echo ${accession}
ls ${accession}
if [[ ${bait_completed::1} == "T" ]]; then
  accession=${accession%?????????}
  echo ${accession}
  mkdir -p results
  mkdir -p results/01-bait-mhc
  mv ./${accession}.fastq.gz results/01-bait-mhc/${accession}.fastq.gz
else
  accession=${accession%????}
  echo ${accession}
  tar -xvf ${accession}.tar
fi
ls results/01-bait-mhc/${accession}.fastq.gz
tar -zxf ref.tar.gz
rm ./._*
snakemake --snakefile ./iwes_fast.smk \
--cores ${cores} \
--configfile ./config.yaml \
--config accession=${accession}
rm -f ./*_R1_001.fastq.gz
rm -f ./*_R2_001.fastq.gz
snakemake --snakefile ./iwes_fl.smk \
--cores ${cores} \
--configfile ./config.yaml \
--config accession=${accession}

rm ./config.yaml
rm ./ref.tar.gz
rm ./iwes_fl.smk
rm ./iwes_fast.smk

# these get merged later and are redundant
rm -rf ./results/04-converted
# these were created so snakemake has a "loop" to go through and are redundant
rm -rf ./results/07-ipd_ref-fl
######
# move files to the front so they are easier to find.
#######
# this file may want to be deleted in the future as it is 100mb, or resorted to compress it down
# for now it is good to use to investigate why something was missed
mv ./results/10-merged-fl/${accession}.merged.bam ./
# This is a core file
mv ./results/12-chimera-fl/${accession}.filtered.merged.bam ./
mv ./results/13-genotypes-fl/${accession}.genotypes.csv ./
mv ./results/05-merged/${accession}.merged.bam ./diag.${accession}.merged.bam
mv ./results/06-depth/${accession}.depth.txt.gz ./diag.${accession}.depth.txt.gz
mv ./results/11-depth-fl/${accession}.depth.txt.gz ./${accession}.depth.txt.gz
mv ./results/01-bait-mhc/${accession}.read_count.txt ./${accession}.read_count.txt
# zip the rest of the results files so the automated program will return the back up results.
tar -zcvf ${accession}.results.tar.gz results