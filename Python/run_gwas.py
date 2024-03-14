# Program to run GWAS for Ath-HMDP metabolome phenotypes, including 4304 untargeted metabolites
# USAGE: run_gwas_metabolome.py -i <input dir> -o <output dir> 
# prints to STDOUT

import os
import glob
import pylab
import datetime
import argparse
import subprocess

import numpy as np
import pandas as pd
import fastlmm.util.util as flutil

from fastlmm.association import single_snp

def run_gwas(pheno):
    """
    main function to run GWAS, given the input directory containing individual phenotype file
    """
    
    print(pheno)
    print(datetime.datetime.now())
    
    pheno_fn = pheno+'.pheno.txt'
    bed_fn = pheno
    cov_fn = pheno+'.cov.txt'
    
    # run GWAS
    results_df = single_snp(test_snps = bed_fn, pheno = pheno_fn, covar=cov_fn, count_A1=False)
    
    # save result
    results_df.to_csv("{}.gwas".format(pheno), sep="\t")

def main(args):

    run_gwas(args.pheno)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-p', '--pheno',
                        help='phenotype name',
                        type=str,
                        default='')

    args = parser.parse_args()

    main(args)
