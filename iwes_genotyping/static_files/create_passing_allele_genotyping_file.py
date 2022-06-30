# #!/usr/bin/env python
# from argparse import ArgumentParser, ArgumentTypeError
# import pandas as pd
# import pysam
#
# def str2bool(v):
#     if isinstance(v, bool):
#         return v
#     if v.lower() in ('yes', 'true', 't', 'y', '1'):
#         return True
#     elif v.lower() in ('no', 'false', 'f', 'n', '0'):
#         return False
#     else:
#         raise ArgumentTypeError('Boolean value expected.')
#
#
# def parse_process_arrays_args(parser: ArgumentParser):
#     """Parses the python script arguments from bash and makes sure files/inputs are valid"""
#     parser.add_argument('--depth_input_path',
#                         type=str,
#                         help='relative path to depth of coverage file (zipped)',
#                         required=True)
#     parser.add_argument('--merged_bam_path',
#                         type=str,
#                         help='relative path to the merged bam file',
#                         required=True)
#     parser.add_argument('--filtered_allele_list_outpath',
#                         type=str,
#                         help='relative path to the filtered_allele_list',
#                         required=True)
#     parser.add_argument('--filtered_merged_bam_outpath',
#                         type=str,
#                         help='relative path to the filtered_merged_bam_',
#                         required=True)
#     parser.add_argument('--depth_threshold',
#                         type=int,
#                         default=10,
#                         help='depth of coverage requirement for the center of the read',
#                         required=False)
#     parser.add_argument('--edge_distance_threshold',
#                         type=int,
#                         help='distance from edge to consider depth of coverage in positions',
#                         default=0,
#                         required=False)
#     parser.add_argument('--maximum_start_position_gap',
#                         type=int,
#                         help='how far apart the start position can be between reads for the chimera criteria',
#                         default=30,
#                         required=False)
#     parser.add_argument('--minimum_bridge_read_length',
#                         type=int,
#                         help='How mapped (untrimmed) portion of hte read must be to be considered',
#                         default=70,
#                         required=False)
#
#
# def get_process_arrays_args():
#     """	Inputs arguments from bash
#     Gets the arguments, checks requirements, returns a dictionary of arguments
#     Return: args - Arguments as a dictionary
#     """
#     parser = ArgumentParser()
#     parse_process_arrays_args(parser)
#     args = parser.parse_args()
#     return args
#
#
# args = get_process_arrays_args()
# # same arguments to a local variable by same name as the argument
# passing_diag_allele_path = args.passing_diag_allele_path
# passing_fl_allele_path = args.passing_fl_allele_path
# genotypes_outpath = args.genotypes_outpath
# filtered_merged_bam_outpath = args.filtered_merged_bam_outpath
#
# input:
# "results/12-chimera-fl/" + accession + ".allele_list.csv",
# "results/06-depth/" + accession + ".allele_list.csv"
#
# output:
# "results/13-genotypes-fl/" + accession + ".genotypes.csv",
#
# import pandas as pd
# df = pd.read_csv(input[0], header=None, names=['allele'])
# df_diag =pd.read_csv(input[1], header=None, names=['allele'])
# print(df_diag)
# diag_allele_list = list(df_diag['allele'].unique())
# ipd_allele_list = list(df['allele'].unique())
# missing_diag_list = []
# for diag_i in diag_allele_list:
#     should_present_allele_list = ipd_diag_dict[diag_i]
#     at_least_one_allele = any(item in should_present_allele_list for item in ipd_allele_list)
#     if not at_least_one_allele:
#         missing_diag_list.append(diag_i)
# df['read_ct'] = 1
# df['accession'] = accession
# if len(missing_diag_list) > 0 :
#     df_diag_missing = pd.DataFrame({'allele':missing_diag_list})
#     df_diag_missing['read_ct'] = 0.99
#     df_diag_missing['accession'] = accession
#     df = pd.concat([df, df_diag_missing], ignore_index=True)
# df.to_csv(output[0],index=False)