bait_completed=F
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

#if [[ $bait_completed == "T" ]]; then
#  mkdir -p results
#  mkdir -p results/01-bait-mhc
#  mv ./${accession}.fastq.gz results/01-bait-mhc/${accession}.fastq.gz
#fi
#

# tar -zxvf ref.tar.gz
#tar -xvf ${accession}.tar
snakemake --snakefile /Users/dabaker3/github/iwes_genotyping/static_files/iwes_fast.smk \
--cores ${cores} \
--configfile ./config.yaml \
--config accession=${accession}

snakemake --snakefile /Users/dabaker3/github/iwes_genotyping/static_files/iwes_fl.smk \
--cores ${cores} \
--configfile ./config.yaml \
--config accession=${accession}

# tar -zcvf {accession}.results.tar.gz results