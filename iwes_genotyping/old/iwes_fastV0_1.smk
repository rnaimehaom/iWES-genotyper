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
accession = config['accession']
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
    accession = config["accession"]
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
        touch("results/" + accession + "_diag_finished.txt")

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
        shell("bbmap.sh in={input[0]} int=t outm={output[0]} ref={input[1]} semiperfectmode=t threads=1 nodisk=t >/dev/null 2>&1")

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
    ZIP compress depth file
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
        minimum_bridge_read_legth=70
    threads: 1
    run:
        import pandas as pd
        import pysam
        # open the depth file
        df = pd.read_csv(input[0],header=None, sep='\t', names=['allele', 'position', 'depth'])
        # figure out the
        df['MAX_POSITION'] = df.groupby(['allele'])['position'].transform('max')
        df["DISTANCE_FROM_EDGE"] = [y - x if y/2 < x else x - 1 for x, y in zip(df['position'],df['MAX_POSITION'])]
        df['DEPTH_THRESHOLD'] = [0 if (x <= params.depth_threshold) & (y > params.edge_distance_threshold) else 1 for
                                 x, y in zip(df['depth'],df['DISTANCE_FROM_EDGE'])]
        # diag_edge_distance_threshold = 0
        # diag_depth_threshold = 10
        df_summary = df[['allele', 'DEPTH_THRESHOLD']].groupby(['allele'])['DEPTH_THRESHOLD'].min().reset_index()
        df_summary_pass = df_summary[df_summary['DEPTH_THRESHOLD'] > 0]
        allele_list = list(df_summary_pass['allele'].unique())
        passing_alleles=[]
        infile = pysam.AlignmentFile(input[1], "rb")
        bam_header = infile.header.copy().to_dict()
        primer_group_unique_list = []
        outfile_groups = {}
        i = 1

        ref_dict = {}
        ref_name_old = 'ZXAZXV'
        position_old = 0
        compare_start_old = 0
        for segment in infile:
            if i % 100000 == 0:
                print(i)
            i += 1

            compare_start = segment.reference_start
            compare_end = segment.reference_end
            if (compare_end is None) or (compare_start is None) :
                continue
            # cigartuples = segment.cigartuples
            ref_name = segment.reference_name
            if ref_name in allele_list:
                if ref_name == ref_name_old:
                    if compare_start - compare_start_old > params.maximum_start_position_gap:
                        print(compare_start_old,compare_start )
                        ref_dict[ref_name] = False

                else:
                    ref_dict[ref_name] = True
                if compare_end - compare_start > params.minimum_bridge_read_legth:
                    compare_start_old = compare_start
                ref_name_old = ref_name

        infile.close()

        for key, value in ref_dict.items():
            if value:
                passing_alleles.append(key)
        df_passing = pd.DataFrame({'allele':passing_alleles})
        df_passing.to_csv(output[0], header=False, index=False)
        infile = pysam.AlignmentFile(input[1], "rb")
        bam_header = infile.header.copy().to_dict()

        outfile = pysam.AlignmentFile(output[1], "wb", header=bam_header)
        i = 1
        for segment in infile:
            if i % 100000 == 0:
                print(i)
            i += 1
            ref_name = segment.reference_name
            if ref_name in passing_alleles:
                outfile.write(segment)
        infile.close()
        outfile.close()

rule create_allele_list_fl:
    """
    exhautively map baited reads by mapping one at a time
    to reference sequences that have exons 2-4
    using bbmap
    suppress output stderr and stdout because it consumes a lot of unnecessary space
    """
    input:
        "results/01-bait-mhc/" + accession + ".fastq.gz",
        "results/06-depth/" + accession + ".allele_list.csv",
    output:
        'results/07-allele-list/' + accession + '.allele_list_fl.txt',
        'results/07-allele-list/' + accession + '.allele_list_fl_num.txt'
    threads: 1
    run:
        exhaustive_result = 'results/02-exhaustive-map-fl'
        sorted_bam_result = 'results/03-sorted-bam-fl'
        sample=os.path.basename(input[1])[:-6]
        df_missing = pd.read_csv(input[2], sep='\t', header=None, names=['allele'])
        missing_allele_list = list(df_missing['allele'].unique())
        # ipd_num_dict
        included =False
        df = pd.read_csv(input[1], sep='\t', header=None, names=['allele'])
        diag_present_list = list(df['allele'].unique())
        # print(diag_present_list)
        ipd_allele_list = []
        # print(diag_present_list)
        for diag_i in diag_present_list:
            if diag_i in ipd_diag_dict.keys():
                ipd_allele_list = ipd_allele_list + ipd_diag_dict[diag_i]
        # print(ipd_allele_list)
        os.makedirs(exhaustive_result, exist_ok=True)
        bam_list = []
        ipd_allele_list = ipd_allele_list + missing_allele_list
        with open(output[0], 'w') as f:
            f.write('\n'.join(ipd_allele_list))
        ipd_num_list = [ipd_num_dict[x] for x in ipd_allele_list]
        with open(output[1], 'w') as f:
            f.write('\n'.join(ipd_num_list))
