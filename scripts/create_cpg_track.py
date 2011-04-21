#!/usr/bin/env python

"""
create_cpg_track.py

Create BED track for CpG sites. The track file can be uploaded to the
UCSC Genome Browser.
"""

import sys
import argparse
import os
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import HUMAN_CHRS, MOUSE_CHRS
from MethylAnalyzer.UtilityFuncs import check_file

def main(chrs, cpg_dir, out_dir):
    LOGFH = open('create_cpg_tracks.LOG', 'w')
    LOGFH.write('Start the program ...  %s\n' % datetime.datetime.now())
    # Get track values for each chromosome
    for chr in chrs:
        LOGFH.write('... Read CpGs on [ %s ] ...  %s\n' % (chr, datetime.datetime.now()))
        # Some fixed info in the output file
        header1 = 'browser position ' + chr + ':1-10000'
        header_cpg = 'track name="'+chr+' CpG" description="CpGs on '+chr+'" color=0,0,0'
        header_re = 'track name="'+chr+' RE" description="RE sites on '+chr+'" color=255,0,0'
        header_mcrbc = 'track name="'+chr+' McrBC" description="McrBC sites on '+chr+'" color=0,0,255'
        # Read CpG file
        try:
            cpgfile = check_file(chr+'_cpgs', cpg_dir)
        except MethError, e:
            print >> stderr, e.value
            sys.exit(2)
        infh = open(cpgfile)
        cpgout = open(os.path.join(out_dir, chr+'_cpgs.bed'), 'w')
        reout = open(os.path.join(out_dir, chr+'_re.bed'), 'w')
        mcrbcout = open(os.path.join(out_dir, chr+'_mcrbc.bed'), 'w')
        cpgout.write(header1 + '\n')
        cpgout.write(header_cpg + '\n')
        reout.write(header1 + '\n')
        reout.write(header_re + '\n')
        mcrbcout.write(header1 + '\n')
        mcrbcout.write(header_mcrbc + '\n')
        for line in infh:
            line_list = line.rstrip().split('\t')
            position = int(line_list[1])
            isre = int(line_list[2])
            ismcrbc = int(line_list[3])
            cpgout.write('%s\t%d\t%d\n' % (chr, position, position+2))
            if isre:
                reout.write('%s\t%d\t%d\n' % (chr, position, position+2))
            if ismcrbc:
                mcrbcout.write('%s\t%d\t%d\n' % (chr, position, position+2))
        infh.close()
        cpgout.close()
        reout.close()
        mcrbcout.close()
    LOGFH.write('Finish the program ...  %s\n' % datetime.datetime.now())

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create BED track for CpG sites of the human or mouse genome. \
    The track file can be uploaded to the UCSC Genome Browser.')
    parser.add_argument('--genome', help='human or mouse genome, default=human', default='human')
    parser.add_argument('--outdir', help='directory storing BED track files for CpG sites, default=current dir', default=os.getcwd())
    parser.add_argument('cpg_dir', help='directory storing CpG annotation files created by parse_sites.py')
    # Parse arguments
    args = parser.parse_args()
    if args.genome == 'human':
        chrs = HUMAN_CHRS
    elif args.genome == 'mouse':
        chrs = MOUSE_CHRS
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid genome, human/mouse genome available'
        sys.exit(2)
    if os.path.isdir(args.cpg_dir) is True and os.path.isdir(args.outdir) is True:
        cpg_dir = args.cpg_dir
        out_dir = args.outdir
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid directories'
        sys.exit(2)
    main(chrs, cpg_dir, out_dir)
    
