#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')
def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--out_dir',
                        type=str,
                        help='directory to output the results',
                        required=True)
    parser.add_argument('--allele_haplotype_caller_path',
                        type=str,
                        help='premade allele haplotype caller path',
                        required=True)
    parser.add_argument('--animal_lookup_path',
                        type=str,
                        help='animal lookup path from labkey',
                        required=True)
    parser.add_argument('--in_dir',
                        type=str,
                        help='where the *.genotypes.csv files are',
                        required=True)
    parser.add_argument('--submit_name',
                        type=str,
                        help='name of the submit file',
                        required=True)
def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
out_dir = args.out_dir
allele_haplotype_caller_path = args.allele_haplotype_caller_path
animal_lookup_path = args.animal_lookup_path
in_dir = args.in_dir
submit_name = args.submit_name

import pandas as pd
import os
import numpy as np


# out_dir = '/Volumes/T7/MHC_genotyper/baylor_31_ipd_diag'
# allele_haplotype_caller_path = '/Volumes/T7/MHC_genotyper/ref/ipd_disagnostic_wes_mamu_allele_haplotype.csv'
# animal_lookup_path = '/Volumes/T7/MHC_genotyper/baylor_31/animal_lookup_baylor31.csv'
# in_dir = '/Volumes/T7/MHC_genotyper/baylor_31_test/results_2/08-genotypes-fl'
# submit_name = 'submit_name'
xlsx_filepath = None
os.makedirs(out_dir, exist_ok=True)


def error_code(x, y, z):
    if y > 2:
        return 'TMH:({0})'.format(z)
    if y == 2:
        return x
    if y == 1:
        # add a dash for a second option, as it needs to be looked at to some degree
        return [x, '-']
    return "NO HAPLO"


# this is a bit complicated:
# First, we need unique names for columns that will/will not be merged
df_ref_lookup = pd.read_csv(allele_haplotype_caller_path)
df_ref_lookup.rename(columns={'allele': 'allele_diagnostic'}, inplace=True)
# Next we want the count of alleles, if it has this count it meets the "all" criteria
df_counts = df_ref_lookup[['HAPLOTYPE_CALL', 'TYPE', 'allele_diagnostic']].drop_duplicates()
df_counts = df_counts.groupby(['HAPLOTYPE_CALL', 'TYPE'])['allele_diagnostic'].count().reset_index()
df_counts.rename(columns={'allele_diagnostic': 'diag_count'}, inplace=True)
df_ref_lookup = df_ref_lookup.merge(df_counts, on=['HAPLOTYPE_CALL', 'TYPE'])


genotyping_pathlist = os.listdir(in_dir)
genotyping_pathlist = [os.path.join(in_dir, x) for x in genotyping_pathlist if
                       (not x.startswith('._') and x.endswith('.genotypes.csv'))]

# need to untar the files first fo the genotyps
#
df = pd.DataFrame()
for path_i in genotyping_pathlist:
    df_i = pd.read_csv(path_i)
    df = pd.concat([df, df_i], ignore_index=True)
df.rename(columns={'accession': 'gs_id'}, inplace=True)
df['allele_idp'] = ['-'.join(x.split('-')[1:]) for x in df['allele']]
df['MHC_TYPE'] = [x.split('_')[0] if y == 1 else x.split('_')[1] for x, y in zip(df['allele'], df['read_ct'])]
type_dict = {'Mamu-A1': 'MHC_A_HAPLOTYPES', 'Mamu-A2': 'MHC_A_HAPLOTYPES', 'Mamu-B': 'MHC_B_HAPLOTYPES'}
#  'Mamu-A3':'MHC_A_HAPLOTYPES', #'Mamu-A4', 'Mamu-A6', 'Mamu-AG1',
# 'Mamu-AG2', 'Mamu-AG3', 'Mamu-AG4', 'Mamu-AG5', 'Mamu-AG6',
# 'Mamu-B02Ps', 'Mamu-B10Ps', 'Mamu-B11L', 'Mamu-B14Ps', 'Mamu-B17',
# 'Mamu-B21Ps',
# 'Mamu-E', 'Mamu-G', 'Mamu-I', 'Mamu-J'}

df['TYPE'] = [type_dict[x] if x in type_dict.keys() else x for x in df['MHC_TYPE']]
df_diag = df[df['read_ct'] == .99]
df_ref_diag = df_ref_lookup[['HAPLOTYPE_CALL', 'allele_diagnostic', 'diag_count', 'TYPE']].drop_duplicates()
allele_idp_list = []
allele_diag_long_list = []
for idx, row in df_diag.iterrows():

    allele_temp_list = [x for x in df_ref_diag['allele_diagnostic'] if (x in row['allele'])]
    allele_idp_list.append(allele_temp_list)
    # allele_temp_list = [row['allele_diagnostic'] for x in df_diag['allele_diagnostic'] if (row['allele_diagnostic'] in x)]
    if len(allele_temp_list) > 0:
        allele_diag_long_list.append(allele_temp_list[0])
    else:
        allele_diag_long_list.append('')
df_diag['match'] = [len(x) for x in allele_idp_list]
df_diag['allele_diagnostic'] = allele_diag_long_list
df_diag_f = df_diag[df_diag['match'] > 0]
df_diag_f = df_diag_f.merge(df_ref_diag, on=['TYPE', 'allele_diagnostic'])
df_diag_f.drop(columns=['match'], inplace=True)

df_merge = df.merge(df_ref_lookup[['TYPE', 'allele_idp', 'HAPLOTYPE_CALL', 'allele_diagnostic', 'diag_count']],
                    on=['TYPE', 'allele_idp'], how='inner')
print(df_diag_f)

df_merge = pd.concat([df_merge, df_diag_f], ignore_index=True)

df_i = df_merge.copy()
df_counts_i = df_i[['gs_id', 'allele_diagnostic', 'TYPE', 'HAPLOTYPE_CALL']].drop_duplicates()
df_counts_i = df_counts_i.groupby(['gs_id', 'HAPLOTYPE_CALL', 'TYPE'])['allele_diagnostic'].count().reset_index()
df_counts_i.rename(columns={'allele_diagnostic': 'sample_counts'}, inplace=True)

df_merge_dedup = df_merge[['gs_id', 'HAPLOTYPE_CALL', 'TYPE', 'diag_count']].drop_duplicates()
df_temp = df_merge_dedup.merge(df_counts_i, on=['gs_id', 'HAPLOTYPE_CALL', 'TYPE'])
df_haplo_pass = df_temp[df_temp['diag_count'] <= df_temp['sample_counts']]
df_haplo_pass['PASS_COUNT'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['HAPLOTYPE_CALL'].transform('count')

df_haplo_pass['rank'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['diag_count'].transform('rank', method='first')

df_haplo_pass['TMH'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['HAPLOTYPE_CALL'].transform('; '.join)

df_haplo_pass['HAPLOTYPE_FINAL'] = [error_code(x, y, z) for x, y, z in
                                    zip(df_haplo_pass['HAPLOTYPE_CALL'], df_haplo_pass['PASS_COUNT'],
                                        df_haplo_pass['TMH'])]
# if it has only one haplotype, it will make a list [haplotype_call, '-'] which then can be xploded
df_haplo_pass = df_haplo_pass.explode('HAPLOTYPE_FINAL')
df_haplo_pass['rank'] = df_haplo_pass.groupby(['gs_id', 'TYPE'])['diag_count'].transform('rank', method='first')
df_haplo_pass['TYPE'] = ['{0} {1}'.format(x.replace('_', ' '), int(y)) for x, y in
                         zip(df_haplo_pass['TYPE'], df_haplo_pass['rank'])]
df_haplo_pass = df_haplo_pass[df_haplo_pass['rank'] < 3]
# gs_id_cols = list(df_haplo_pass['gs_id'].unique())
df_haplo_pivot = df_haplo_pass.pivot_table(values='HAPLOTYPE_FINAL',
                                           index=['TYPE'],
                                           columns=['gs_id'],
                                           aggfunc=np.max,
                                           fill_value='NO HAPLO').reset_index()

df_haplo_pivot.rename_axis(['index'], inplace=True, axis=1)
print(df_haplo_pivot)
## get the path name of xlsx (or create one)
## then concat the genotype data to the haplotype calls and export to excel

if xlsx_filepath is None:
    xlsx_filepath = os.path.join(out_dir, '{0}.pivot.xlsx'.format(submit_name))
df['# Obs'] = df.groupby('allele_idp')['read_ct'].transform('count')
genotype_pivot = df.pivot_table(values='read_ct',
                                index=['allele_idp', '# Obs'],
                                columns=['gs_id'],
                                aggfunc=np.max,
                                fill_value='').reset_index()
genotype_pivot.rename_axis(['index'], inplace=True, axis=1)
genotype_pivot.rename(columns={'allele_idp': 'gs_id'}, inplace=True)
df_haplo_pivot.rename(columns={'TYPE': 'gs_id'}, inplace=True)

df_allele_count_summary = df.groupby('gs_id')['read_ct'].count().reset_index()
df_allele_count_summary.rename(columns={'read_ct': '# Alleles Identified'}, inplace=True)
df_count_summary_pivot = df_allele_count_summary.pivot_table(values='# Alleles Identified',
                                                             columns=['gs_id'],
                                                             aggfunc=np.max,
                                                             fill_value=0).reset_index()
df_count_summary_pivot.rename_axis(['index2'], inplace=True, axis=1)
df_count_summary_pivot.rename(columns={'index': 'gs_id'}, inplace=True)
df_count_summary_pivot.rename_axis(['index'], inplace=True, axis=1)


df_gs_id = pd.DataFrame([list(df_haplo_pivot.columns)], columns=list(df_haplo_pivot.columns))
df_xlx_pivot = pd.concat([df_gs_id, df_count_summary_pivot, df_haplo_pivot, genotype_pivot], ignore_index=True)

df_xlx_pivot.rename(columns={'gs_id': 'Animal IDs'}, inplace=True)
if os.path.exists(animal_lookup_path):
    animal_lookup = pd.read_csv(animal_lookup_path)
    animal_lookup_dict = {}
    for idx, row in animal_lookup.iterrows():
        animal_lookup_dict[row['gs_id']] = row['animal_id']
    gs_id_list = list(df_xlx_pivot.columns)
    df_xlx_pivot.rename(columns=animal_lookup_dict, inplace=True)
    animal_list = list(animal_lookup['animal_id'])
    animal_list2 = []

    for gs_id_i in gs_id_list:
        if gs_id_i in animal_lookup_dict.keys():
            animal_list2.append(animal_lookup_dict[int(gs_id_i)])
    animal_list2.sort()
    col_order = ['Animal IDs', '# Obs'] + animal_list2
else:
    col_names = list(df_xlx_pivot.columns)
    col_names.remove('Animal IDs')
    col_names.remove('# Obs')
    col_order = ['Animal IDs', '# Obs'] + col_names
df_xlx_pivot = df_xlx_pivot[col_order]
df_xlx_pivot.to_excel(xlsx_filepath, sheet_name='Sheetname_1', index=False)

