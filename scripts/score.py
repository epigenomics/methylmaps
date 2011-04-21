#!/usr/bin/env python

"""
score.py

Estimate DNA methylation of CpGs
"""

import sys
import argparse
import os
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import MethStateParser
from MethylAnalyzer.UtilityFuncs import check_chr, check_file

def main(p_bar, chr, chr_len, meth_data, out_dir):
    LOGFH = open(os.path.join(out_dir, 'score.LOG'), 'a')    
    LOGFH.write('Start the program for [%s] ...  %s\n' % (os.path.basename(meth_data), str(datetime.datetime.now())))
    LOGFH.flush()
    # 1. Rebuild CpG sites
    LOGFH.write('... Rebuild CpGs ...  %s\n' % str(datetime.datetime.now()))
    LOGFH.flush()
    methparser = MethStateParser(chr_len, 1000)
    methparser.rebuild_cpgs(meth_data)
    # 2. Estimate methylation probability and write to the output
    LOGFH.write('... Estimate methylation probability ...  %s\n' % str(datetime.datetime.now()))
    try:
        p1_p2, p_bar_new = methparser.estimate_methprob(p_bar)
    except MethError, e:
        LOGFH.write('MethError: [ %s ] \n' % e.value)
        LOGFH.write('No methylation probabilities estimated ...\n\n')
        LOGFH.close()
        sys.exit(2)
    LOGFH.write("... p1/p2 = %.2f, p_bar_new = %.2f ...\n" % (p1_p2, p_bar_new))
    LOGFH.flush()
    cpg_sites = methparser.get_cpgs()
    # 3. Save methylation probabilities
    outfile = os.path.join(out_dir, 'methstatus_' + chr)
    outfh = open(outfile, 'w')
    for cpg in cpg_sites.sorted_iter():
        record = cpg.get_tablerecord(chr)
        outfh.write('%s\t%.2f\t%.2f\n' % (record, cpg.p, cpg.n))
    outfh.close()
    LOGFH.write('Finish the program ...  %s\n\n' % str(datetime.datetime.now()))
    LOGFH.close()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Estimate DNA methylation states for CpGs in a sample')
    parser.add_argument('--out_dir', help='directory for output files, default=currect dir', default=os.getcwd())
    parser.add_argument('meth_ave', type=float, help='global methylation level estimated by LUMA')
    parser.add_argument('chr', help='chromosome')
    parser.add_argument('chr_len', type=int, help='chromosome length')
    parser.add_argument('methdata', help='methylation data generated by filter.py')
    # Parse arguments
    args = parser.parse_args()
    p_bar = args.meth_ave
    chr_len = args.chr_len
    try:
        chr = check_chr(args.chr)
        meth_data = check_file(args.methdata)
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, e.value
        sys.exit(2)
    if os.path.isdir(args.out_dir) is True:
        out_dir = os.path.abspath(args.out_dir)
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid directories'
        sys.exit(2)
    main(p_bar, chr, chr_len, meth_data, out_dir)
