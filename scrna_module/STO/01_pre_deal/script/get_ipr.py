#!/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3
import sys
import os
import h5py
import argparse

__doc__ = '''
处理ipr命名
'''
__author__ = 'leiguo'
__mail__ = 'leiguo@genome.cn'

def ipr(ipr, tar, outdir):
    tar_name = os.path.basename(tar).replace('.tar.gz', '')
    f = h5py.File(ipr)
    qc_result_file = f["ssDNA/ImageInfo"].attrs['QCResultFile']
    if tar_name != qc_result_file:
        tar_name = qc_result_file
    cmd = 'cp -rf {} {}/{}.ipr'.format(ipr, outdir, tar_name)
    cmd += ' && cp -rf {} {}/{}.tar.gz'.format(tar, outdir, tar_name)
    os.system(cmd)
    f.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}'.format(__author__, __mail__))
    parser.add_argument('-i', '--ipr', help='QC result ipr', dest='ipr', required=True)
    parser.add_argument('-t', '--tar', help='QC result tar', dest='tar', required=True)
    parser.add_argument('-o', '--outdir', help='outdir', dest='outdir', required=True)
    args = parser.parse_args()

    ipr(args.ipr, args.tar, args.outdir)


if __name__ == '__main__':
    main()

