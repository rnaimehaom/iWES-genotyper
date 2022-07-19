
for arg in "$@"; do
  shift
  case "$arg" in
    "--gdrive_out_dir") set -- "$@" "-g" ;;
    "--in_dir") set -- "$@" "-i" ;;
    "--out_dir") set -- "$@" "-o" ;;
    "--submit_name") set -- "$@" "-s" ;;
    "--allele_haplotype_caller_path") set -- "$@" "-h" ;;
    "--animal_lookup_path") set -- "$@" "-a" ;;
    "--cassowary_path") set -- "$@" "-c" ;;
    *) set -- "$@" "$arg"
  esac
done

while getopts g:o:i:h:a:s:c: opt ; do
  case $opt in
    g) gdrive_out_dir=$OPTARG ;;
    i) in_dir=$OPTARG ;;
    o) out_dir=$OPTARG ;;
    h) allele_haplotype_caller_path=$OPTARG ;;
    a) animal_lookup_path=$OPTARG ;;
    s) submit_name=$OPTARG ;;
    c) cassowary_path=$OPTARG ;;
    *) usage; exit 1;;
  esac
done

cd ${in_dir}
LATEST_DIR=`find ${in_dir} -maxdepth 1 -type d -name "iwes*" -print0 | xargs -0 ls -d | tail -n 1`
if [[ "${LATEST_DIR}" == "" ]]; then
    echo "no iwes dir"
    exit
fi
echo "${LATEST_DIR}"

# PANGOLIN_REPORT_DIR=${LATEST_DIR}/pangolin_reports
mkdir -p ${LATEST_DIR}
cd ${in_dir}
find ${LATEST_DIR} -maxdepth 2 -type f -name "*.tar.gz" > genotypes_list.txt
cd ${LATEST_DIR}
for tar_file in `cat ../genotypes_list.txt`; do
  echo ${tar_file}
  tar -xzvf ${tar_file} '*.genotypes.csv'
  tar -xzvf ${tar_file} '*.read_count.txt'
done
find ./ -name '*.genotypes.csv' -exec mv '{}' ./ \;
find ./ -name '*.read_count.txt' -exec mv '{}' ./ \;
script_directory=$(dirname "$0")
echo ${script_directory}
echo "python3 ${script_directory}/aggregate_haplotype.py --submit_name=${submit_name} \
      --out_dir=${out_dir} \
      --allele_haplotype_caller_path=${allele_haplotype_caller_path} \
      --animal_lookup_path=${animal_lookup_path} \
      --in_dir=${LATEST_DIR}
"

python3 ${script_directory}/aggregate_haplotype.py --submit_name=${submit_name} \
--out_dir=${out_dir} \
--allele_haplotype_caller_path=${allele_haplotype_caller_path} \
--animal_lookup_path=${animal_lookup_path} \
--in_dir=${LATEST_DIR}

echo "current dir ${out_dir}"
echo "coping files to ${gdrive_out_dir}/${submit_name}"
echo "coping files to ${cassowary_path}/${submit_name}"

mkdir -p "${gdrive_out_dir}/${submit_name}"
mkdir -p "${cassowary_path}/${submit_name}"
rsync -aP ${out_dir}/ "${gdrive_out_dir}/${submit_name}/"
rsync -aP ${out_dir}/ "${cassowary_path}/${submit_name}/"
rsync -aP ${LATEST_DIR}/*.tar.gz "${gdrive_out_dir}/${submit_name}/"
rsync -aP ${LATEST_DIR}/*.tar.gz "${cassowary_path}/${submit_name}/"

python3 ${script_directory}/notify_slack.py ${cassowary_path}/${submit_name}/ ${submit_name}


