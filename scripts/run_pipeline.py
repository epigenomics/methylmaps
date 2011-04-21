#!/usr/bin/env python

"""
run_pipeline.py

The master script to generate all required sub-scripts for the analysis pipeline.
"""

import sys
import os
import re
import argparse
import datetime

from MethylAnalyzer.UtilityFuncs import check_file
from MethylAnalyzer.MethError import MethError

# Pipeline scripts
PIPELINE = {1: 'parse_mates.py', 2: 'filter.py', 3: 'score.py'}

def main(parafile, refile, mcrbcfile, out_dir, RUN):
    ## Get parameters required for the pipeline
    para_dict = get_paras(parafile)
    chrlen_dict = get_lens(para_dict['CHR_LENGTH'])
    ## Step 1. Parse *.mates files
    # 1.1 Make fragment directory
    fragdir = os.path.join(out_dir, 'fragments')
    os.mkdir(fragdir)
    # 1.2 Make scripts
    step1_re = ' '.join([PIPELINE[1], '--out_dir', fragdir, 're', para_dict['CMAP'], refile])
    step1_mcrbc = ' '.join([PIPELINE[1], '--out_dir', fragdir, 'mcrbc', para_dict['CMAP'], mcrbcfile])
    # 1.3 Run scripts or save them to the script file
    if RUN:
        print '**** Start step 1 - parse *.mates files ****', datetime.datetime.now()
        os.system(step1_re)
        os.system(step1_mcrbc)
        print '**** End step 1 ****\n', datetime.datetime.now()
    else:
        script_file = os.path.join(fragdir, 'scripts_step1')
        step1fh = open(script_file, 'w')
        step1fh.write('%s\n%s\n' % (step1_re, step1_mcrbc))
        step1fh.close()
    ## Step 2. Filter RE/McrBC fragments
    # 2.1 Make filtering directory
    fltdir = os.path.join(out_dir, 'filtering')
    os.mkdir(fltdir)
    # 2.2 Make scripts
    step2_scripts = []
    for chr in chrlen_dict.keys():
        refrag = os.path.join(fragdir, chr+'_re.bed')
        mcrbcfrag = os.path.join(fragdir, chr+'_mcrbc.bed')
        script = ' '.join([PIPELINE[2], '--out_dir', fltdir, chr, para_dict['ANNO_DIR'], refrag, mcrbcfrag])
        step2_scripts.append(script)
    # 2.3 Run scripts or save them to the script file
    if RUN:
        print '**** Start step 2 - filter RE/McrBC fragments ****', datetime.datetime.now()
        for script in step2_scripts:
            os.system(script)
        print '**** End step 2 ****\n', datetime.datetime.now()
    else:
        script_file = os.path.join(fltdir, 'scripts_step2')
        step2fh = open(script_file, 'w')
        for script in step2_scripts:
            step2fh.write('%s\n' % script)
        step2fh.close()
    ## Step 3. Estimate DNA methylation probabilities
    # 3.1 Make scoring directory
    scoredir = os.path.join(out_dir, 'methscores')
    os.mkdir(scoredir)
    # 3.2 Make scripts
    step3_scripts = []
    for chr in chrlen_dict.keys():
        methfile = os.path.join(fltdir, 'methdata_'+chr)
        script = ' '.join([PIPELINE[3], '--out_dir', scoredir, para_dict['METH_AVE'], chr, chrlen_dict[chr], methfile])
        step3_scripts.append(script)
    # 3.3 Run scripts or save them to the script file
    if RUN:
        print '**** Start step 3 - estimate methylation probabilities ****'
        for script in step3_scripts:
            os.system(script)
        print '**** End step 3 ****\n'
    else:
        script_file = os.path.join(scoredir, 'scripts_step3')
        step3fh = open(script_file, 'w')
        for script in step3_scripts:
            step3fh.write('%s\n' % script)
        step3fh.close()
            
        
def get_lens(lenfile):
    """
    Get chromosome lengths
    File format
       chr   length
    """
    chrlen_dict = {}
    lines = open(lenfile).readlines()
    for line in lines:
        chr, length = line.rstrip().split('\t')
        chrlen_dict[chr] = length
    return chrlen_dict

def get_paras(parafile):
    """
    Get methylation parameters for the analysis pipeline.
    Format of the parameter file (separated by TAB):
       CMAP        /Users/user/...
       ANNO_DIR    /Users/user/...
       METH_AVE    /Users/user/...
       CHR_LENGTH  /Users/user/...
    """
    para_dict = {}
    lines = open(parafile).readlines()
    for line in lines:
        if not re.search('^#', line):  # skip comments
            k, v = line.rstrip().split('\t')
            para_dict[k] = v
    return para_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The master script to generate all sub-scripts for the data analysis pipeline')
    parser.add_argument('para', help='parameters required in the pipeline')
    parser.add_argument('re_mates', help='mates file for RE fragments')
    parser.add_argument('mcrbc_mates', help='mates file for McrBC fragments')
    parser.add_argument('out_dir', help='directory for all output files of the pipeline')
    parser.add_argument('--run', type=bool, nargs='?', const=True, default=False, \
                        help='run the analysis pipeline or just save scripts, default=False')
    # Parse arguments
    args = parser.parse_args()
    try:
        parafile = check_file(args.para)
        refile = check_file(args.re_mates)
        mcrbcfile = check_file(args.mcrbc_mates)
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, e.value
        sys.exit(2)
    if os.path.isdir(args.out_dir):
        out_dir = os.path.abspath(args.out_dir)
    main(parafile, refile, mcrbcfile, out_dir, args.run)
    
