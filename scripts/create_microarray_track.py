#!/usr/bin/env python

"""
create_microarray_track.py

Create methylation profile tracks in the microarray format for the UCSC Genome Browser
"""

import sys
import argparse
import os
import datetime

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.MethBuilder import SiteParser, HUMAN_MAIN_CHRS, MOUSE_MAIN_CHRS
from MethylAnalyzer.UtilityFuncs import check_file

def main(chrs, scale, step, meth_dir, cpg_dir, out_dir, sname):
    LOGFH = open('create_microarray_tracks.LOG', 'w')
    LOGFH.write('Start the program ...  %s\n' % datetime.datetime.now())
    # Create methylation tracks for each chromosome
    for chr in chrs:
        # 1. Build CpGs
        LOGFH.write('... Build CpG sites for [ %s ] ...  %s\n' % (chr, datetime.datetime.now()))
        try:
            cpgfile = check_file(chr+'_cpgs', PATH=cpg_dir)
        except MethError, e:
            print >> sys.stderr, e.value
            continue
        siteparser = SiteParser(chr)
        siteparser.parse_sites(cpgfile, 'CpG')
        cpg_sites = siteparser.get_sites()
        for cpg in cpg_sites:  # define empty meth score
            cpg.meth_score = '-10000'
        # 2. Get methylation scores
        # note files are saved in 'methstatus_chrN'
        LOGFH.write('... Obtain methylation scores ...\n')
        try:
            mfile = check_file('methstatus_' + chr, PATH=meth_dir)
        except MethError, e:
            print >> sys.stderr, e.value
            sys.exit(2)
        mfh = open(mfile)
        for line in mfh:
            line_list = line.rstrip().split('\t')
            coordinate = int(line_list[1])
            p = float(line_list[8]) * 2 * scale - scale  # convert to a new score based on scales
            # Assign methylation score to the CpG
            cpg_sites[coordinate].meth_score = str(round(p, 2))
        mfh.close()
        # 3. Write to the output file
        LOGFH.write('... Save the track file ...\n')
        # Some fixed info in the output file
        header = 'track type="array" expScale=' + str(scale) + ' expStep=' + str('%.1f' % step) + ' expNames="' + \
                 sname + '," name="' + chr + '" description="Microarray custom track for ' + chr + '"' 
        fix1 = '\t'.join(['meth', '500', '+'])
        fix2 = '\t'.join(['0', '1', '1,', '0,', '1', '0,'])
        outfile = 'meth_' + chr + '.bed'
        fout = open(outfile, 'w')
        fout.write(header + '\n')
        for coord in cpg_sites.sorted_coords():
            cpg = cpg_sites[coord]
            record = '\t'.join([chr, str(cpg.coordinate), str(cpg.coordinate+1), fix1, str(cpg.coordinate),
                                str(cpg.coordinate+1), fix2, cpg.meth_score])
            fout.write(record + '\n')
        fout.close()        
        del siteparser
        del cpg_sites
    LOGFH.write('Finsh the program ...  %s\n' % datetime.datetime.now())
    LOGFH.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create custom track files in the microarray format based on \
    methylation scores.')
    parser.add_argument('--genome', help='human or mouse genome, default=human', default='human')
    parser.add_argument('--scale', help='scale used in the header of the track file, default=2', default=2)
    parser.add_argument('--step', help='color key step value in the header, default=0.2', default=0.2)
    parser.add_argument('--out_dir', help='directory for created custom track files, default=current dir', default=os.getcwd())
    parser.add_argument('sample', help='sample name')
    parser.add_argument('meth_dir', help='directory for methylation score files')
    parser.add_argument('cpg_dir', help='directory for CpG annotation files')
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
    if os.path.isdir(args.meth_dir) is True and os.path.isdir(args.cpg_dir) is True and \
           os.path.isdir(args.out_dir) is True:
        meth_dir = args.meth_dir
        cpg_dir = args.cpg_dir
        out_dir = args.out_dir
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid directories'
        sys.exit(2)
    main(chrs, args.scale, args.step, meth_dir, cpg_dir, out_dir, args.sample)
