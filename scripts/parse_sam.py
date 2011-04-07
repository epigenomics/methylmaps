#!/usr/bin/env python

"""
Parse SAM/BAM file to find proper read pairs
"""

####
# If the SAM/BAM is corrected created by a specific sequencing platform,
#   the easiest way to extract paired reads is to check the flags.
#   if 0x0002, then the read is mapped in a proper pair.
####
# Specific criteria to find paired read in a SAM/BAM file
# 1. same chromosome
# 2. same strand
# 3. acceptable insert size (0-15Kb)
# 4. coorect read ordering
#    4.1 for plus strand: R3 pos < F3 pos
#    4.2 for minus strand: F3 pos < R3 pos
####
# Relative positions of reads
#   '+' strand:  ----R3----     ----F3----
#                |        |     |        |
#              .pos    .aend  .pos     .aend
#
#   '-' strand:  ----F3----     ----R3----
#                |        |     |        |
#              .pos    .aend  .pos     .aend
#                |                       |
#           left coord              right coord
####

import os
import sys
import pysam
import pdb
import argparse
import datetime
import bisect
from string import lower

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../methylanalyzer'))
from MethError import MethError
from UtilityFuncs import check_file

def main(FLAG, fformat, readlib, chr, infile, out_dir):
    LOGFH = open(os.path.join(out_dir, 'parse_sam.LOG'), 'a')
    LOGFH.write('Start to parse [%s] reads in the SAM/BAM file for [%s] ...  %s\n' % (readlib, infile, str(datetime.datetime.now())))
    LOGFH.flush()
    # Read the input file
    samfile = pysam.Samfile(infile, fformat)
    chrs = samfile.references  # chromosomes
    if chr not in chrs:
        print >> sys.stderr, 'Wrong chromosome name'
        sys.exit(2)
    # Parse and write to the output files organized by chromosome
    LOGFH.write('... Parse properly paired reads ...  %s\n' % str(datetime.datetime.now()))
    LOGFH.flush()
    outfh = generate_fh(chr, readlib, out_dir)
    ids_sorted = []
    for read in samfile:
        if chrs[read.tid] != chr:
            print >> sys.stderr, 'Contain reads on [%s] instead of [%s]' % (chrs[read.tid], chr)
            sys.exit(2)
        # Check if the read has been saved to the output
        insert_point = bisect.bisect_right(ids_sorted, read.qname)
        if ids_sorted[insert_point-1:insert_point] == [read.qname]:
            continue
        # Check if the read is paired
        if FLAG:
            is_paired = ispaired_byflag(bin(read.flag))
        else:
            is_paired = ispaired_bydef(read)
        # Conclusion: The read is properly paired and has not been saved
        if is_paired:
            ids_sorted.insert(insert_point, read.qname)                        
            record = get_bedrecord(chr, read)
            outfh.write('%s\n' % record)
            if len(ids_sorted) % 100000 == 0:
                LOGFH.write('...... [%d] properly paired reads parsed ......\n' % len(ids_sorted))
                LOGFH.flush()
    samfile.close()
    outfh.close()
    LOGFH.write('Finish the program ...  %s\n\n' % (str(datetime.datetime.now())))
    LOGFH.close()

def ispaired_byflag(flag_bin):
    "The second field of flag records if the read is mapped in a proper pair"
    if flag_bin[-2] == '1':
        return True
    else:
        return False

def ispaired_bydef(read):
    "Use defined criteria to check if two reads are properly paired"
    # Step 1. on same chromosome
    if read.tid != read.mrnm:
        return False
    # Step 2: on same strand
    if read.is_reverse != read.mate_is_reverse:
        return False
    # Steps 3 & 4. acceptable insert size and correct read ordering
    if read.is_read1:
        f3_pos = read.pos
        r3_pos = read.mpos
    else:
        f3_pos = read.mpos
        r3_pos = read.pos
    frag_len = abs(f3_pos - r3_pos)
    if (not read.is_reverse) and (not read.mate_is_reverse) and r3_pos < f3_pos and frag_len > 0 and frag_len < 15000:  # on '+' strand
        return True
    elif read.is_reverse and read.mate_is_reverse and f3_pos < r3_pos and frag_len > 0 and frag_len < 15000:  # on '-' strand
        return True
    else:
        return False

def get_bedrecord(chr, read):
    """
    Save fragment information based on F3 and R3 reads
    """
    if read.is_read1:
        f3_pos = read.pos
        f3_len = read.rlen
        r3_pos = read.mpos
        r3_len = read.qlen
    else:
        f3_pos = read.mpos
        f3_len = read.qlen
        r3_pos = read.pos
        r3_len = read.rlen
    if read.is_reverse:
        strand = '-'
        left_coord = f3_pos  # 0-based
        right_coord = r3_pos + r3_len  # 1-based
    else:
        strand = '+'
        left_coord = r3_pos  # 0-based
        right_coord = f3_pos + f3_len  # 1-based
    score = '1'
    itemRgb = '0'
    blockCount = '2'
    blockSize = '25,25,'
    blockStarts = '0,' + str(right_coord - left_coord - 25)
    mates_id = read.qname
    record = '\t'.join([chr, str(left_coord), str(right_coord), mates_id, score,
                        strand, str(left_coord), str(right_coord), itemRgb,
                        blockCount, blockSize, blockStarts])
    return record

def generate_fh(chr, readlib, out_dir):
    "Generate the output file handle"
    fname = os.path.join(out_dir, chr + '_' + lower(readlib) + '.bed')
    fh = open(fname, 'w')
    if readlib == 'RE':  # red
        color = '255,0,0,'
    else:  # McrBC, blue
        color = '0,0,255,'
    # Print headers
    fh.write('browser position %s:1-10000\n' % chr)
    fh.write('track name="%s %s" description="Fragment data from the %s library" color=%s\n' % \
            (chr, readlib, readlib, color))
    return fh

def get_readlib(arg):
    lib_dict = {'re': 'RE', 'mcrbc': 'McrBC'}
    return lib_dict[arg]

def get_format(arg):
    'Get file open mode based on input file format'
    fdict = {'sam': 'r', 'bam': 'rb'}
    return fdict[arg]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse mate-pair reads for a specific chromosome. Reads are saved in SAM/BAM format.',
                                     epilog='Save coordinates of fragments that are formed by properly paired reads \
                                     in BED formate, generating files like chr1_re.bed.')
    parser.add_argument('--flag', type=bool, nargs='?', const=True, default=False, \
                        help='use flag value to extract paired reads, default=False')
    parser.add_argument('--out_dir', default=os.getcwd(), help='directory for parsed fragments, default=current dir')
    parser.add_argument('format', choices=['sam', 'bam'], help='input file format')
    parser.add_argument('lib', choices=['re', 'mcrbc'], help='sequences generated by the RE or McrBC library')
    parser.add_argument('chr', help='chromosome name, e.g., chr1')
    parser.add_argument('input', help='input file name')
    # Parse arguments
    args = parser.parse_args()
    readlib = get_readlib(args.lib)
    fformat = get_format(args.format)
    try:
        infile = check_file(args.input)
    except MethError, e:
        parser.print_usage()
        print >> sys.stderr, 'MethError: ', e.value
        sys.exit(2)
    if os.path.isdir(args.out_dir):
        out_dir = os.path.abspath(args.out_dir)
    else:
        parser.print_usage()
        print >> sys.stderr, 'Invalid output directory'
        sys.exit(2)
    main(args.flag, fformat, readlib, args.chr, infile, out_dir)
        
