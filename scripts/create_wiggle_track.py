#!/usr/bin/env python

"""
create_wiggle_track.py

Create methylation profile tracks in the wiggle format for the UCSC Genome Browser
"""

import sys
import os
import argparse
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import SiteParser, HUMAN_MAIN_CHRS, MOUSE_MAIN_CHRS
from MethylAnalyzer.UtilityFuncs import check_file

def main(chrs, meth_dir, out_dir, sname):
    LOGFH = open('create_wiggle_tracks.LOG', 'w')
    LOGFH.write('Start the program ...  %s\n' % datetime.datetime.now())
    # Create methylation tracks for each chromosome
    for chr in chrs:
        LOGFH.write('... Obtain methylation scores for [ %s ] ...\n' % chr)
        # 1. Get methylation scores
        # note files are saved in 'methstatus_chrN'
        try:
            mfile = check_file('methstatus_' + chr, PATH=meth_dir)
        except MethError, e:
            print >> sys.stderr, e.value
            sys.exit(2)
        cpg_sites = {}  # coordinate: score
        mfh = open(mfile)
        for line in mfh:
            line_list = line.rstrip().split('\t')
            coordinate = int(line_list[1])
            cpg_sites[coordinate] = line_list[8]
        mfh.close()
        # 2. Write to the output file
        # Some fixed info in the output file
        header1 = 'track type=wiggle_0 name=' + chr + ' description=Wiggle custom track for ' \
                 + sname + '_' + chr + ' color=128,0,0 visibility=full'
        header2 = 'variableStep chrom=' + chr
        outfile = 'meth_' + chr + '.wig'
        fout = open(outfile, 'w')
        fout.write(header1 + '\n' + header2 + '\n')
        sorted_coords = cpg_sites.keys()
        sorted_coords.sort()
        for coord in sorted_coords:
            record = '\t'.join([str(coord), cpg_sites[coord]])
            fout.write(record + '\n')
        fout.close()        
        del cpg_sites
    LOGFH.write('Finsh the program ...  %s\n' % datetime.datetime.now())
    LOGFH.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create custom tracks in the \
    wiggle format based on methylation scores.')
    parser.add_argument('--genome', help='human or mouse genome, default=human', \
                        default='human')
    parser.add_argument('--out_dir', help='directory for created custom track files, default=current dir', default=os.getcwd()) 
    parser.add_argument('sample', help='sample name')
    parser.add_argument('meth_dir', help='directory for methylation score files')
    # Parse arguments
    args = parser.parse_args()
    if args.genome == 'human':
        chrs = HUMAN_MAIN_CHRS
    elif args.genome == 'mouse':
        chrs = MOUSE_MAIN_CHRS
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid genome, human/mouse genome available'
        sys.exit(2)
    if os.path.isdir(args.meth_dir) is True and os.path.isdir(args.out_dir) is True:
        meth_dir = args.meth_dir
        out_dir = args.out_dir
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid directories'
        sys.exit(2)
    main(chrs, meth_dir, out_dir, args.sample)

