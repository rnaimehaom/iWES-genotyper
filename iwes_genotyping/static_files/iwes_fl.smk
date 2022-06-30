from pathlib import Path
from Bio import SeqIO
import os
import pandas as pd
import json
import shutil
## CONFIG ##

# configfile: "config/config.yaml"


## INSPECT INPUT FASTQ ##

# test if FASTQ files exist in `fastq_folder`
# if so, use these files
# otherwise, download FASTQ from SRA
# get FASTQ folder and extension from config
fastq_folder = config["fastq"]["fastq_folder"]
fastq_extension = config["fastq"]["R1_fastq_extension"]
accession = str(config['accession'])
fastq_source = ""
if any(File.endswith(fastq_extension) for File in os.listdir(fastq_folder)):
    # use FASTQ files that are already downloaded
    # set accession to basename of FASTQ file
    for path in Path(fastq_folder).glob("*" + fastq_extension):
        accession = str(path).split("_")[0]
    fastq_source = "local"

    # paths to use for baiting with local FASTQ
    R1_fastq = fastq_folder + "/" + accession + "_R1_001.fastq.gz"
    R2_fastq = fastq_folder + "/" + accession + "_R2_001.fastq.gz"
else:
    # specify SRA accession in snakemake invokation
    accession = str(config['accession'])
    R1_fastq = ""
    R2_fastq = ""
    fastq_source = "sra"

missing_alleles = config["missing_alleles"]
mapping_reference_fasta = config["mapping_reference_fasta"]
diag_ref_fasta = config["diag_ref_fasta"]
indvl_allele_fl_fa_dir = config["indvl_allele_fl_fa_dir"]
indvl_allele_fa_dir = config["indvl_allele_fa_dir"]
ref_fasta, = glob_wildcards(indvl_allele_fa_dir + "/{allele}.fasta")
# clean ._ files
ref_fasta = [x for x in ref_fasta if not x.startswith('._')]
concatenated_fasta_name = os.path.basename(diag_ref_fasta)
filepath = 'results/07-allele-list/' + accession + '.allele_list_fl_num.txt'
with open('ref/diag_to_ipd_lookup.json') as f_in:
    ipd_diag_dict = json.load(f_in)
with open('ref/ipd_num_lookup.json') as f_in:
    ipd_num_dict = json.load(f_in)

with open(filepath) as f:
    # ) and (x in ipd_diag_dict.keys())
    allele_list = f.read().split('\n')

per_sample_list = [x.strip('\"')  for x in allele_list if len(x) > 0]
ipd_ref_fl = 'results/07-ipd_ref-fl'
os.makedirs(ipd_ref_fl, exist_ok=True)
for sample_i in per_sample_list:
    fp_dest = os.path.join(ipd_ref_fl,'{0}.fasta'.format(sample_i))
    fp_src = os.path.join(indvl_allele_fl_fa_dir,'{0}.fasta'.format(sample_i))
    shutil.copyfile(fp_src, fp_dest)
# clean ._ files
wildcard_constraints:
    allele_fl = "(?!\._).*"

rule all:
    input:
        "results/12-chimera-fl/" + accession + ".filtered.merged.bam",
        "results/13-genotypes-fl/" + accession + ".genotypes.csv"
    run:
        shell("touch results/{accession}_finished.txt")

rule exhaustive_mapping_fl:
    """
    exhautively map baited reads by mapping one at a time
    to reference sequences that have exons 2-4
    using bbmap
    suppress output stderr and stdout because it consumes a lot of unnecessary space
    """
    input:
        "results/01-bait-mhc/" + accession + ".fastq.gz",
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta',
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta.fai'
    output:
        temp('results/07-exaustive-fl/' + '{allele_fl}.sam')
    threads: 1
    run:
        shell("bbmap.sh -Xmx4g in={input[0]} int=t outm={output[0]} ref={input[1]} semiperfectmode=t threads=1 nodisk=t >/dev/null 2>&1")

rule sort_sam_fl:
    """
    sort SAM files following bbmap
    make list of output files to use with samtools merge
    """
    input:
        'results/07-exaustive-fl/' + '{allele_fl}.sam'
    output:
        temp('results/08-sortsam-fl/' + '{allele_fl}.sorted.sam')
    threads: 1
    run:
        shell("samtools sort {input[0]} -o {output[0]}")


rule index_fasta_fl:
    """
    create index for each individual FASTA file
    sometimes fails for no reason, so retry three times
    """
    input:
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta'
    output:
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta.fai'
    threads: 1
    run:
        shell("samtools faidx {input[0]}")

rule convert_bam_fl:
    """
    convert SAM file to BAM
    """
    input:
        'results/08-sortsam-fl/' + '{allele_fl}.sorted.sam',
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta',
        'results/07-ipd_ref-fl/' + '{allele_fl}.fasta.fai'
    output:
        temp('results/09-cv-bam-fl/' + '{allele_fl}.sorted.bam')
    run:
        shell("samtools view -b -h -T {input[1]} -o {output[0]} {input[0]} \
        && samtools index {output[0]}")

rule merge_bam_fl:
    """
    merge sorted SAM files into single SAM file
    """
    input:
        expand("results/09-cv-bam-fl/{per_sample}.sorted.bam",
            per_sample=per_sample_list)
    output:
        "results/10-merged-fl/" + accession + ".merged.bam"
    params:
        per_sample_files=lambda wildcards, input: " ".join(input)
    run:
        shell("find ./results/09-cv-bam-fl  -name '*.bam' > ./results/10_filelist.txt")
        shell("mkdir -p ./results/10_bam_split")
        shell("split -l 200 ./results/10_filelist.txt ./results/10_bam_split/bam_segment")
        shell("find ./results/10_bam_split -name 'bam_segment*' > ./results/10_split_files.txt")
        shell("for f in `cat ./results/10_split_files.txt`; do \
        samtools merge ${{f}}.bam -b ${{f}} && samtools index ${{f}}.bam; \
        done")
        shell("find ./results/10_bam_split -name '*.bam' > ./results/10_merged_bam_files.txt")
        shell("samtools merge {output[0]} -b ./results/10_merged_bam_files.txt && samtools index {output[0]}")
        shell("rm -rf ./results/10_bam_split")

rule save_ref_with_bam_fl:
    """
    save the reference sequence with the BAM file for inspection in Geneious
    """
    input:
        mapping_reference_fasta
    output:
        "results/10-merged-fl/ipd_2022_iwes.fasta"
    run:
        ref_name = os.path.basename(mapping_reference_fasta)
        ref_dest = os.path.join(sorted_bam_result, ref_name)
        shell("cp {input[0]} {output[0]}")


rule compute_depth_fl:
    """
    use samtools depth to compute depth at each position of SAM file
    by default several types of reads are filtered from depth calculations
    add back with -g [UNMAP,SECONDARY,QCFAIL,DUP]
    """
    input:
        "results/10-merged-fl/" + accession + ".merged.bam",
        mapping_reference_fasta
    output:
        temp("results/11-depth-fl/" + accession + ".depth.txt"),
    threads: 1
    run:

        shell("samtools depth -a {input[0]} -g UNMAP,SECONDARY,QCFAIL,DUP -o {output[0]}")

rule compress_depth_fl:
    """
    ZIP compress depth file
    """
    input:
        "results/11-depth-fl/" + accession + ".depth.txt",
    output:
        "results/11-depth-fl/" + accession + ".depth.txt.gz",
    threads: 1
    run:
        shell("gzip -c {input[0]} > {output[0]}")

rule filter_depth_chimeras_fl:
    """
    ZIP compress depth file
    """
    input:
        "results/11-depth-fl/" + accession + ".depth.txt.gz",
        "results/10-merged-fl/" + accession + ".merged.bam"
    output:
        "results/12-chimera-fl/" + accession + ".allele_list.csv",
        "results/12-chimera-fl/" + accession + ".filtered.merged.bam",
        "results/12-chimera-fl/" + accession + '.finished.txt'
    params:
        edge_distance_threshold=20,
        depth_threshold=10,
        maximum_start_position_gap=28,
        minimum_bridge_read_length=70
    threads: 1
    run:
        shell("python3 ./filter_depth_chimera.py --depth_input_path={input[0]} \
--merged_bam_path={input[1]} \
--filtered_allele_list_outpath={output[0]} \
--filtered_merged_bam_outpath={output[1]} \
--edge_distance_threshold={params.edge_distance_threshold} \
--depth_threshold={params.depth_threshold} \
--maximum_start_position_gap={params.maximum_start_position_gap} \
--minimum_bridge_read_length={params.minimum_bridge_read_length}")
        shell('touch {output[2]}')

rule genotype:
    """
    identify SAM mappings where there is complete read support at least specified depth
    """
    input:
        "results/12-chimera-fl/" + accession + ".allele_list.csv",
        "results/06-depth/" + accession + ".allele_list.csv"

    output:
        "results/13-genotypes-fl/" + accession + ".genotypes.csv",

    threads: 1
    run:
        import pandas as pd
        df = pd.read_csv(input[0], header=None, names=['allele'])
        df_diag =pd.read_csv(input[1], header=None, names=['allele'])
        print(df_diag)
        diag_allele_list = list(df_diag['allele'].unique())
        ipd_allele_list = list(df['allele'].unique())
        missing_diag_list = []
        for diag_i in diag_allele_list:
            should_present_allele_list = ipd_diag_dict[diag_i]
            at_least_one_allele = any(item in should_present_allele_list for item in ipd_allele_list)
            if not at_least_one_allele:
                missing_diag_list.append(diag_i)
        df['read_ct'] = 1
        df['accession'] = accession
        if len(missing_diag_list) > 0 :
            df_diag_missing = pd.DataFrame({'allele':missing_diag_list})
            df_diag_missing['read_ct'] = 0.99
            df_diag_missing['accession'] = accession
            df = pd.concat([df, df_diag_missing], ignore_index=True)
        df.to_csv(output[0],index=False)
