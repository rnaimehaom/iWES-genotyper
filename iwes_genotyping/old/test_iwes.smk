# from pathlib import Path
from Bio import SeqIO
# from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from pathlib import Path
# import shutil
# import scripts
# from scripts import prepare_ipd_reference
# from scripts import dho_utils
import json
import os
import pandas as pd
import json

with open('ref/ipd_to_diag_lookup.json') as f_in:
    ipd_diag_dict = json.load(f_in)
# cd /Volumes/T7/MHC_genotyper/baylor_31_test/
# snakemake --snakefile /Users/dabaker3/github/iwes_genotyping/test_iwes.smk --cores 1
individual_allele_fl_fasta_folder = 'ref/2022-02-02_ipd_fl'

ref_fasta_fl, = glob_wildcards(individual_allele_fl_fasta_folder + "/{allele}.fasta")
# ref_fasta_fl='ref/ipd-mhc-mamu-2022-02-02_fl.txt'
missing_alleles='ref/missing_ipd_alleles.csv'
mapping_reference_fasta = 'ref/ipd-mhc-mamu-2022-02-02_fl.fasta'
accession='102875'
rule all:
    input:
        "results/06-depth-fl/" + accession + ".depth.txt.gz"
    run:
        touch("results/"++ accession + "finished.txt")
rule exhaustive_mapping_fl:
    """
    exhautively map baited reads by mapping one at a time
    to reference sequences that have exons 2-4
    using bbmap
    suppress output stderr and stdout because it consumes a lot of unnecessary space
    # filter list of alleles
    # It takes snake make about 5-10 seconds to run an empty file.
    # This is a long time when you have 2000 files.
    # therefore we will not thread this step since there are only ~100 alleles per sample  
    """
    input:
        "results/01-bait-mhc/" + accession + ".fastq.gz",
        "results/06-depth/" + accession + ".allele_list.csv",
        'ref/missing_ipd_alleles.csv'
    output:
        temp("results/02-exhaustive-map-fl/{allele}.sam"),
    threads: 1
    run:
        sample=os.path.basename(input[1])[:-6]
        df_missing = pd.read_csv(input[2], header=None, names=['allele'])
        missing_allele_list = list(df_missing['allele'].unique())
        name = None
        fasta_sequences = SeqIO.parse(open(input[1]), 'fasta')
        for fasta in fasta_sequences:
            name= fasta.id
            break
        included =False
        if name not in missing_allele_list:
            df = pd.read_csv(input[2], header=None, names=['allele'])
            diag_present_list = list(df['allele'].unique())
            if name in ipd_diag_dict.keys():
                diag_list = ipd_diag_dict[name]
                if any(item in diag_list for item in diag_present_list):
                    included = True
        if included or (name in missing_allele_list):
            shell("bbmap.sh in={input[0]} int=t outm={output[0]} ref={input[1]} semiperfectmode=t threads=1 nodisk=t >/dev/null 2>&1")
        shell('touch {output[0]}')
#             shell('echo "@HD	VN:1.4	SO:unsorted\
# @SQ	SN:{sample}	LN:1\
# @PG	ID:BBMap	PN:BBMap	VN:38.24	CL:java -ea -Xmx3200m align2.BBMap build=1" > {output[0]}')

rule sort_sam_fl:
    """
    sort SAM files following bbmap
    make list of output files to use with samtools merge
    """
    input:
        "results/02-exhaustive-map-fl/{allele}.sam",
    output:
        temp("results/03-sorted-fl/{allele}.sorted.sam"),
    threads: 1
    run:
        filesize = os.path.getsize(input[0])
        if filesize != 0:
            shell("samtools sort {input[0]} -o {output[0]}")
        shell('touch {output[0]}')
rule index_fasta_fl:
    """
    create index for each individual FASTA file
    sometimes fails for no reason, so retry three times
    """
    input:
        individual_allele_fl_fasta_folder + "/{allele}.fasta"
    output:
        individual_allele_fl_fasta_folder + "/{allele}.fasta.fai"
    threads: 1
    run:
        filesize = os.path.getsize(input[0])
        if filesize != 0:
            shell("samtools faidx {input[0]}")
        shell('touch {output[0]}')
rule convert_bam_fl:
    """
    convert SAM file to BAM
    """
    input:
        "results/03-sorted-fl/{allele}.sorted.sam",
        individual_allele_fl_fasta_folder + "/{allele}.fasta",
        individual_allele_fl_fasta_folder + "/{allele}.fasta.fai"
    output:
        temp("results/04-converted-fl/{allele}.sorted.bam")
    run:
        filesize = os.path.getsize(input[0])
        if filesize != 0:
            shell("samtools view -b -h -T {input[1]} -o {output[0]} {input[0]} \
            && samtools index {output[0]}")
        shell('touch {output[0]}')

rule merge_bam_fl:
    """
    merge sorted SAM files into single SAM file
    """
    input:
        expand("results/04-converted-fl/{per_sample}.sorted.bam",per_sample=ref_fasta_fl)
    output:
        "results/05-merged-fl/" + accession + ".merged.bam"
    params:
        per_sample_files=lambda wildcards, input: " ".join(input)
    run:
        shell("samtools merge {output[0]} results/04-converted-fl/*.bam && samtools index {output[0]}")

# rule save_ref_with_bam_fl:
#     """
#     save the reference sequence with the BAM file for inspection in Geneious
#     """
#     input:
#         mapping_reference_fasta
#     output:
#         "results/04-converted-fl/" + concatenated_fasta_name
#     run:
#         shell("cp {input[0]} {output[0]}")


rule compute_depth_fl:
    """
    use samtools depth to compute depth at each position of SAM file
    by default several types of reads are filtered from depth calculations
    add back with -g [UNMAP,SECONDARY,QCFAIL,DUP]
    """
    input:
        "results/05-merged-fl/" + accession + ".merged.bam",
        mapping_reference_fasta
    output:
        temp("results/06-depth-fl/" + accession + ".depth.txt"),
    threads: 1
    run:

        shell("samtools depth -a {input[0]} -g UNMAP,SECONDARY,QCFAIL,DUP -o {output[0]}")

rule compress_depth_fl:
    """
    ZIP compress depth file
    """
    input:
        "results/06-depth-fl/" + accession + ".depth.txt",
    output:
        "results/06-depth-fl/" + accession + ".depth.txt.gz",
    threads: 1
    run:
        shell("gzip -c {input[0]} > {output[0]}")
