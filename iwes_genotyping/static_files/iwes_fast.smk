from pathlib import Path
from Bio import SeqIO
import os
import pandas as pd
import json
import re

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
# results/01-bait-mhc/${accession}.fastq.gz
with open('ref/diag_to_ipd_lookup.json') as f_in:
    ipd_diag_dict = json.load(f_in)
with open('ref/ipd_num_lookup.json') as f_in:
    ipd_num_dict = json.load(f_in)

# clean ._ files
wildcard_constraints:
    allele = "(?!\._).*"

rule all:
    input:
        'results/07-allele-list/' + accession + '.allele_list_fl.txt'
    run:

        shell("touch results/{accession}_diag_finished.txt")

if not os.path.exists("results/01-bait-mhc/" + accession + ".fastq.gz"):
    rule bait_mhc:
        """
        baiting extracts MHC reads from iWES FASTQ
        this makes it faster to do exhaustive searching for MHC-mapped reads that map perfectly
        it also simplifies troubleshooting and optimizing genotyping parameters
        download FASTQ from SRA with fasterq-dump and send to stdout, then
        map FASTQ to CY0333 gDNA MHC reads
        save mapped reads as FASTQ file
        """
        input:
            mapping_reference_fasta
        output:
            "results/01-bait-mhc/" + accession + ".fastq.gz"
        params:
            accession=accession,
            extra="--skip-technical",
            r1_fastq=R1_fastq,
            r2_fastq=R2_fastq
        threads: 8
        run:
            if fastq_source == "sra":
                # if fastq_source is SRA
                shell("fasterq-dump {params.accession} \
                --threads {threads} --split-spot --stdout -p  \
                | bbmap.sh -Xmx12g in=stdin.fq int=t outm={output[0]} ref={input[0]} semiperfectmode=t threads={threads}")
            if fastq_source == "local":
                # if fastq_source is local
                shell("bbmap.sh -Xmx12g in={params.r1_fastq} in2={params.r2_fastq} outm={output[0]} ref={input[0]} semiperfectmode=t threads={threads}")
# java -ea -Xmx4500m -Xms4500m -cp /opt/conda/opt/bbmap-38.90-3/current/ align2.BBMap build=1 overwrite=true \
#     fastareadlen=500 \
#                  in=results/01-bait-mhc/102923.fastq.gz int=t \
#     outm=results/02-exhaustive-map/572.sam \
#     ref=ref/2022-02-02_diag/572.fasta \
#     semiperfectmode=t threads=1 nodisk=t
rule exhaustive_mapping:
    """
    exhautively map baited reads by mapping one at a time
    to reference sequences that have exons 2-4
    using bbmap
    suppress output stderr and stdout because it consumes a lot of unnecessary space
    """
    input:
        "results/01-bait-mhc/" + accession + ".fastq.gz",
        indvl_allele_fa_dir + "/{allele}.fasta"
    output:
        temp("results/02-exhaustive-map/{allele}.sam"),
    threads: 1
    run:
        shell("bbmap.sh -Xmx4g in={input[0]} int=t outm={output[0]} ref={input[1]} semiperfectmode=t threads=1 nodisk=t")
        # shell("bbmap.sh in={input[0]} int=t outm={output[0]} ref={input[1]} semiperfectmode=t threads=1 nodisk=t >/dev/null 2>&1")
rule sort_sam:
    """
    sort SAM files following bbmap
    make list of output files to use with samtools merge
    """
    input:
        "results/02-exhaustive-map/{allele}.sam",
    output:
        temp("results/03-sorted/{allele}.sorted.sam"),
    threads: 1
    run:
        shell("samtools sort {input[0]} -o {output[0]}")


rule index_fasta:
    """
    create index for each individual FASTA file
    sometimes fails for no reason, so retry three times
    """
    input:
        indvl_allele_fa_dir + "/{allele}.fasta"
    output:
        indvl_allele_fa_dir + "/{allele}.fasta.fai"
    threads: 1
    run:
        shell("samtools faidx {input[0]}")

rule convert_bam:
    """
    convert SAM file to BAM
    """
    input:
        "results/03-sorted/{allele}.sorted.sam",
        indvl_allele_fa_dir + "/{allele}.fasta",
        indvl_allele_fa_dir + "/{allele}.fasta.fai"
    output:
        temp("results/04-converted/{allele}.sorted.bam")
    run:
        shell("samtools view -b -h -T {input[1]} -o {output[0]} {input[0]} \
        && samtools index {output[0]}")

rule merge_bam:
    """
    merge sorted SAM files into single SAM file
    Some nodes/systems have file limits set.  
    By looping over 200 files at a time we should not meet that file open limit.
    The default file limit to be around 1024, but there is a lot opened in the background 
        (and possibly on the shared underlying node) that can cause this to be reached as low as 500 files.
    """
    input:
        expand("results/04-converted/{per_sample}.sorted.bam",per_sample=ref_fasta)
    output:
        "results/05-merged/" + accession + ".merged.bam"
    params:
        per_sample_files=lambda wildcards, input: " ".join(input)
    run:
        shell("find ./results/04-converted  -name '*.bam' > ./results/04_filelist.txt")
        shell("mkdir -p ./results/04_bam_split")
        shell("split -l 200 ./results/04_filelist.txt ./results/04_bam_split/bam_segment")
        shell("find ./results/04_bam_split -name 'bam_segment*' > ./results/04_split_files.txt")
        shell("for f in `cat ./results/04_split_files.txt`; do \
        samtools merge ${{f}}.bam -b ${{f}} && samtools index ${{f}}.bam; \
        done")
        shell("find ./results/04_bam_split -name '*.bam' > ./results/04_merged_bam_files.txt")
        shell("samtools merge {output[0]} -b ./results/04_merged_bam_files.txt && samtools index {output[0]}")
        shell("rm -rf ./results/04_bam_split")

rule save_ref_with_bam:
    """
    save the reference sequence with the BAM file for inspection in Geneious
    """
    input:
        diag_ref_fasta
    output:
        "results/04-converted/" + concatenated_fasta_name
    run:
        shell("cp {input[0]} {output[0]}")

rule compute_depth:
    """
    use samtools depth to compute depth at each position of SAM file
    by default several types of reads are filtered from depth calculations
    add back with -g [UNMAP,SECONDARY,QCFAIL,DUP]
    """
    input:
        "results/05-merged/" + accession + ".merged.bam",
        diag_ref_fasta
    output:
        temp("results/06-depth/" + accession + ".depth.txt"),
    threads: 1
    run:
        shell("samtools depth -a {input[0]} -g UNMAP,SECONDARY,QCFAIL,DUP -o {output[0]}")

rule compress_depth:
    """
    ZIP compress depth file
    """
    input:
        "results/06-depth/" + accession + ".depth.txt",
    output:
        "results/06-depth/" + accession + ".depth.txt.gz",
    threads: 1
    run:
        shell("gzip -c {input[0]} > {output[0]}")

rule filter_depth_chimeras:
    """
    1.) filters the depth of coverage (based on a depth and position from the edge)
    2.) next step filters out chimeras:
    2a.)It ensures there is at least one mapping read a set gap (default 30 positions) apart across the allele
    3.) Next a list of passing alleles is created as a single column csv file
    4.) Next a bam files with mapped reads to the the "passing" alleles 
    """
    input:
        "results/06-depth/" + accession + ".depth.txt.gz",
        "results/05-merged/" + accession + ".merged.bam",
    output:
        "results/06-depth/" + accession + ".allele_list.csv",
        "results/05-merged/" + accession + ".filtered.merged.bam",
    params:
        edge_distance_threshold=0,
        depth_threshold=10,
        maximum_start_position_gap=30,
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

rule create_allele_list_fl:
    """
    exhautively map baited reads by mapping one at a time
    to reference sequences that have exons 2-4
    using bbmap
    suppress output stderr and stdout because it consumes a lot of unnecessary space
    """
    input:
        "results/06-depth/" + accession + ".allele_list.csv",
        missing_alleles
    output:
        'results/07-allele-list/' + accession + '.allele_list_fl.txt',
        'results/07-allele-list/' + accession + '.allele_list_fl_num.txt'
    threads: 1
    run:
        exhaustive_result = 'results/02-exhaustive-map-fl'
        sorted_bam_result = 'results/03-sorted-bam-fl'
        sample=accession
        # open the list of missing alleles do to no corresponding diag region to the fl.
        df_missing = pd.read_csv(input[1], sep='\t', header=None, names=['allele'])
        missing_allele_list = list(df_missing['allele'].unique())
        # ipd_num_dict
        included =False
        # open the list of alleles that past the depth/computational chimera filter
        df = pd.read_csv(input[0], sep='\t', header=None, names=['allele'])
        # get a list of the unique alleles
        diag_present_list = list(df['allele'].unique())
        # print(diag_present_list)
        ipd_allele_list = []
        # print(diag_present_list)
        # convert the list to the corresponding fl sequences, as many of them multi-map.
        for diag_i in diag_present_list:
            if diag_i in ipd_diag_dict.keys():
                ipd_allele_list = ipd_allele_list + ipd_diag_dict[diag_i]
        # print(ipd_allele_list)
        os.makedirs(exhaustive_result, exist_ok=True)
        # add on the missing allele list
        ipd_allele_list = ipd_allele_list + missing_allele_list
        # write teh two lists of alleles and corresponding number list for the fasta files.
        with open(output[0], 'w') as f:
            f.write('\n'.join(ipd_allele_list))
        ipd_num_list = [ipd_num_dict[x] for x in ipd_allele_list]
        with open(output[1], 'w') as f:
            f.write('\n'.join(ipd_num_list))
