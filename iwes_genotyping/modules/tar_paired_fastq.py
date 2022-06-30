#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import subprocess
import os


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
    parser.add_argument('--in_dir',
                        type=str,
                        help='input_directory_to_tar',
                        required=True)
    parser.add_argument('--out_dir',
                        type=str,
                        help='out directory_to_store tar',
                        required=True)
    parser.add_argument('--extension1',
                        type=str,
                        help='forward_extension',
                        default='_R1_001.fastq.gz',
                        required=False)
    parser.add_argument('--extension2',
                        type=str,
                        help='reverse extension',
                        default='_R2_001.fastq.gz',
                        required=False)


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
in_dir = args.in_dir
out_dir = args.out_dir
extension1 = args.extension1
extension2 = args.extension2

# get the list of files
file_list = os.listdir(in_dir)
# get the files with the extension and not the '._' prefix that causes crashes
ext1_list = [x for x in file_list if x.endswith(extension1) and (not x.startswith('._'))]
# convert to the second extension

ext2_list = ['{0}{1}'.format(x[:-len(extension1)], extension2) for x in ext1_list]
# loop through and tar the files in pairs
for ext1_i, ext2_i in zip(ext1_list, ext2_list):
    print(ext1_i, ext2_i)
    # extract the accession number
    accession = ext1_i[:-len(extension1)]
    subprocess.call('tar -cvf {0}/{1}.tar {2} {3}'.format(out_dir, accession, ext1_i, ext2_i),
                    shell=True,
                    cwd=in_dir
                    )
print('Tar in pairs complete')
