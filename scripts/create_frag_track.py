#!/usr/bin/env python

"""
create_frag_track.py

Create custom track files for RE/McrBC fragments (BED format) to
display in the UCSC Genome Browser.
"""

import sys
import argparse
import bisect

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.UtilityFuncs import check_file

def main(header, listfile, fragfile, outfile):
    # Get the id list, doesn't need to be unique
    idlist = []  # [id_length, ...]
    listfh = open(listfile)
    for line in listfh:
        line_list = line.rstrip().split('\t')
        length = str(int(line_list[2]) - int(line_list[1]))
        idlist.append(line_list[0]+'-'+length)
    listfh.close()
    # Sort idlist
    idlist.sort()
    newfh = open(outfile, "w")
    fragfh = open(fragfile)
    # Header line
    if header:
        for i in range(0, header):
            line = fragfh.readline()
            newfh.write(line)
    # The remaining lines
    for line in fragfh:  # BED format
        line_list = line.rstrip().split("\t")
        length = str(int(line_list[2]) - int(line_list[1]))
        id = line_list[3] + '-' + length
        id_insert_point = bisect.bisect_right(idlist, id)
        if id_insert_point != 0 and id_insert_point <= len(idlist) and idlist[id_insert_point-1] == id:
            newfh.write(line)
    fragfh.close()
    newfh.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create custom track files for RE/McrBC fragments (BED format) to \
    display in the UCSC Genome Browser.')
    parser.add_argument('--header', help='number of header lines, default=2', type=int, default=2)
    parser.add_argument('id_list', help='list of [fragment IDs, start coordinate, and end coordinate], \
    generated by filter.py')
    parser.add_argument('frags', help='the origianl fragments (BED format) generated by parse_mates.py')
    parser.add_argument('newfrags', help='the new fragments based on id_list')
    # Parse arguments
    args = parser.parse_args()
    header = args.header
    outfile = args.newfrags
    try:
        listfile = check_file(args.id_list)
        fragfile = check_file(args.frags)
    except MethError, e:
        parser.print_usage()
        print >> stderr, e.value
        sys.exit(2)
    main(header, listfile, fragfile, outfile)
