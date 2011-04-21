#!/usr/bin/env python

"""
parse_sites.py

Parse CpG, RE, and McrBC sites of a given DNA sequence.
"""

import sys
import os
import argparse
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import SiteParser
from MethylAnalyzer.UtilityFuncs import check_chr, check_file

def main(chr, seqfile, outname):
    # Initiate the SiteParser
    LOGFH = open('parse_sites.LOG', 'w')
    LOGFH.write('Start the program for [%s] ...  %s\n' % (chr, str(datetime.datetime.now())))
    LOGFH.flush()
    siteparser = SiteParser(chr)
    # Parse sites
    LOGFH.write('... Start to parse sites ...  %s\n' % str(datetime.datetime.now()))
    LOGFH.flush()
    siteparser.parse_seq(seqfile)
    # Write to the output files
    LOGFH.write('... Save annotation files ...  %s\n' % str(datetime.datetime.now()))
    LOGFH.flush()
    # 1. CpGs
    write2output(outname+'_cpgs', siteparser.sites, chr)
    # 2. RE
    write2output(outname+'_re', siteparser.re_sites, chr)
    # 3. McrBC
    write2output(outname+'_mcrbc', siteparser.mcrbc_sites, chr)
    LOGFH.write('Finish the program for [%s] ...  %s\n' % (chr, str(datetime.datetime.now())))
    LOGFH.close()

def write2output(filename, sites, chr):
    fh = open(filename, 'w')
    for pos in sites.sorted_iter():
        record = pos.get_record(chr)
        fh.write(record + '\n')
    fh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse genomic sequences for CpG, RE, and McrBC sites')
    parser.add_argument('chr', help='chromosome')
    parser.add_argument('fasta', help='Fasta sequence for the chromosome')
    parser.add_argument('--outname', help='name for the output files, default=chrN')
    # Parse arguments
    args = parser.parse_args()
    try:
        chr = check_chr(args.chr)
        seqfile = check_file(args.fasta)
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, e.value
        sys.exit(2)
    if args.outname is not None:
        outname = args.outname
    else:
        outname = chr
    main(chr, seqfile, outname)
                                     
