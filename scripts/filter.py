#!/usr/bin/env python

"""
filter.py

Filter methylated/unmethylated fragments for enzyme recognition sites.
"""

import sys
import argparse
import textwrap
import os
import re
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import SiteParser, FragmentParser, REFilter, MCRBCFilter
from MethylAnalyzer.UtilityFuncs import check_chr, check_file, count_file_lines

DIR = os.path.dirname(__import__('MethylAnalyzer').__file__)

def main(parafile, chr, anno_dir, re_fragfile, mcrbc_fragfile, out_dir):
    # Get annotation files: e.g., chr1_cpgs, chr1_re, chr1_mcrbc
    try:
        cpg_file = check_file(chr+'_cpgs', PATH=anno_dir)
        re_sitefile = check_file(chr+'_re', PATH=anno_dir)
        mcrbc_sitefile = check_file(chr+'_mcrbc', PATH=anno_dir)
    except MethError, e:
        print >> sys.stderr, e.value
        sys.exit(2)
    # 1. Read parameters
    LOGFH = open(os.path.join(out_dir, 'filter.LOG'), 'a')
    LOGFH.write('Start the program for [%s] ...  %s\n' % (chr, str(datetime.datetime.now())))
    LOGFH.flush()
    paras = parse_paras(parafile)
    # 2. Build CpGs
    LOGFH.write('... Build CpG sites ...  %s\n' % str(datetime.datetime.now()))
    siteparser = SiteParser(chr)
    siteparser.parse_sites(cpg_file, 'CpGFull')
    cpg_sites = siteparser.get_sites()
    # 3. Build RE sites first
    LOGFH.write('... Build RE sites ...  %s\n' % str(datetime.datetime.now()))
    siteparser = SiteParser(chr)
    siteparser.parse_sites(re_sitefile, 'RE')
    re_sites = siteparser.get_sites()
    del siteparser
    # 4. Build RE fragments
    LOGFH.write('... Build RE fragments ...  %s\n' % str(datetime.datetime.now()))
    fragparser = FragmentParser()
    # Split fragfile into subsets to avoid memory overflow
    re_linenum = count_file_lines(re_fragfile)
    refh = open(re_fragfile)
    re_fraglines = []
    pre_count = 0
    for count, line in enumerate(refh):
        re_fraglines.append(line)        
        if (count % 10000 == 0 and count != 0) or count + 1 == re_linenum:
            re_frags = [fragment for fragment in fragparser.parse_bedfrags(re_fraglines, 'RE')]
            # 5. Filter RE fragments, cpg_sites are updated
            LOGFH.write('...... Filter RE fragments from lines [%d, %d] ......\n' % (pre_count, count))
            refilter = REFilter(re_frags, cpg_sites, re_sites, paras['outlen'], paras['inlen'])
            refilter.scan()
            # Save information for the passed and failed RE fragments
            write_logfiles(out_dir, re_frags, paras, 're')
            # Reset re_fraglines            
            re_fraglines = []
            del re_frags, refilter
            pre_count = count + 1
    refh.close()
    # 6. Delete RE sites and fragments
    del re_sites, re_fraglines
    # 7. Build McrBC sites
    LOGFH.write('... Build McrBC sites ...  %s\n' % str(datetime.datetime.now()))
    siteparser = SiteParser(chr)
    siteparser.parse_sites(mcrbc_sitefile, 'McrBC')
    mcrbc_sites = siteparser.get_sites()
    del siteparser
    # 8. Build McrBC fragments
    LOGFH.write('... Build McrBC fragments ...  %s\n' % str(datetime.datetime.now()))
    # Split fragfile into subsets to avoid memory overflow
    mcrbc_linenum = count_file_lines(mcrbc_fragfile)
    mcrbcfh = open(mcrbc_fragfile)
    mcrbc_fraglines = []
    pre_count = 0
    for count, line in enumerate(mcrbcfh):
        mcrbc_fraglines.append(line)
        if (count % 10000 == 0 and count != 0) or count + 1 == mcrbc_linenum:
            mcrbc_frags = [fragment for fragment in fragparser.parse_bedfrags(mcrbc_fraglines, 'McrBC')]
            # 9. Filter McrBC fragments, cpg_sites are updated
            LOGFH.write('...... Filter McrBC fragments from lines [%d, %d] ......\n' % (pre_count, count))
            mcrbcfilter = MCRBCFilter(mcrbc_frags, cpg_sites, mcrbc_sites, paras['winlen'], paras['cut_dist'],
                                      paras['min_dist'], paras['max_dist'], paras['mfuzzylen'])
            mcrbcfilter.scan()
            # Save information for the passed and failed McrBC fragments
            write_logfiles(out_dir, mcrbc_frags, paras, 'mcrbc')
            # Reset mcrbc_fraglines
            mcrbc_fraglines = []
            del mcrbc_frags, mcrbcfilter
            pre_count = count + 1
    # 10. Delete McrBC sites and fragments
    del mcrbc_sites, mcrbc_fraglines
    # 11. Save the meth data file
    LOGFH.write('... Write to the output ...  %s\n' % str(datetime.datetime.now()))
    outfile = os.path.join(out_dir, 'methdata_'+chr)
    write_fltresults(outfile, cpg_sites, chr)
    LOGFH.write('Finish the program ...  %s\n\n' % str(datetime.datetime.now()))
    LOGFH.close()

# Methods
def parse_paras(parafile):
    "Parse a set of filtering parameters that are stored in the parameter file"
    # Default values
    paras = {'outlen': 5, 'inlen': 5, 'winlen': 500, 'cut_dist': 50, 'mfuzzylen': 50, 'min_dist': 40,
            'max_dist': 3000, 'failed_2ends': 3, 'failed_mid': 4, 'passed_2ends': 1, 'passed_1end': 2}
    parafh = open(parafile)
    for line in parafh:
        if re.search('^#', line):
            continue
        line_list = line.rstrip().split('\t')
        if line_list[0] in paras.keys():
            paras[line_list[0]] = int(line_list[1])
    parafh.close()
    return paras
               
def write_logfiles(out_dir, frags, paras, lib):
    "Save passed and failed fragments."
    if paras['failed_2ends']:
        fname = os.path.join(out_dir, 'failed_' + lib + 'frags_ends_' + chr)
        write2log(fname, frags, paras['failed_2ends'])
    if paras['failed_mid']:
        fname = os.path.join(out_dir, 'failed_' + lib + 'frags_mid_' + chr)
        write2log(fname, frags, paras['failed_mid'])
    if paras['passed_2ends']:
        fname = os.path.join(out_dir, 'passed_' + lib + 'frags_2ends_' + chr)
        write2log(fname, frags, paras['passed_2ends'])
    if paras['passed_1end']:
        fname = os.path.join(out_dir, 'passed_' + lib + 'frags_1end_' + chr)
        write2log(fname, frags, paras['passed_1end'])
    rfname = os.path.join(out_dir, lib + '_' + chr)
    write2record(rfname, frags, chr)
    
def write2log(filename, frags, filter_code):
    """
    Save fragment ids to a log file for a specific filter_code.
    1 - passed, enzyme site at both ends
    2 - passed, enzyme site at only one end
    3 - failed, enzyme site NOT at either end
    4 - failed, enzyme site NOT in the interior region of the fragment
    """
    fout = open(filename, 'a')
    for fragment in frags:
        if fragment.filter == filter_code:
            record = '\t'.join([fragment.get_matesid(), str(fragment.left_coord), str(fragment.right_coord)])
            fout.write(record + '\n')
    fout.close()

def write2record(filename, frags, chr):
    """
    Save filtering results for all fragments.
    """
    fout = open(filename, 'a')
    for fragment in frags:
        record = fragment.get_filterrecord(chr)
        fout.write(record + '\n')
    fout.close()

def write_fltresults(filename, cpg_sites, chr):
    outfh = open(filename, 'w')
    for cpg in cpg_sites:
        record = cpg.get_methrecord(chr)
        outfh.write(record + '\n')
    outfh.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
                                     Filter methylated/unmethylated fragments for RE or McrBC recognition 
                                     sites at two ends'''),
                                     epilog=textwrap.dedent('''\
                                     The output files for a specific chromosome include
                                         1. failed fragments in four classes: failed_refrags_ends,
                                            failed_refrags_mid, failed_mcrbcfrags_ends, failed_mcrbcfrags_mid
                                         2. passed fragments in four classes: passed_refrags_1end,
                                            passed_refrags_2ends, passed_mcrbcfrags_2ends,
                                            passed_mcrbcfrags_1end
                                         3. save filtering results to files'''))
    parser.add_argument('--para', help='filtering parameters', default=os.path.join(DIR, 'data/filter_para'))
    parser.add_argument('--out_dir', help='directory for output files, default=currect dir', default=os.getcwd())
    parser.add_argument('chr', help='chromosome')
    parser.add_argument('anno_dir', help='directory storing annotation files for CpG, RE, and McrBC sites')
    parser.add_argument('re_frags', help='RE fragments')
    parser.add_argument('mcrbc_frags', help='McrBC fragments')
    # Parse arguments
    args = parser.parse_args()
    try:
        parafile = check_file(args.para)  # parafile
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, e.value
        sys.exit(2)
    if os.path.isdir(args.anno_dir) is True and os.path.isdir(args.out_dir) is True:
        anno_dir = os.path.abspath(args.anno_dir)  # anno_dir
        out_dir = os.path.abspath(args.out_dir)  # out_dir
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid directories'
        sys.exit(2)
    try:
        chr = check_chr(args.chr)  # chr
        re_fragfile = check_file(args.re_frags)  # re_fragfile
        mcrbc_fragfile = check_file(args.mcrbc_frags)  # mcrbc_fragfile
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, e.value
        sys.exit(2)        
    main(parafile, chr, anno_dir, re_fragfile, mcrbc_fragfile, out_dir)
    
