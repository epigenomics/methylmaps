# Utility functions

"""
Utility functions
"""

import os
import re
import fnmatch
from types import IntType, LongType, FloatType, NoneType, StringType
from MethError import MethError

# Define epsilon
EPSILON = 0.0000001

def sortedDictValues(adict):
    "Sort a dict"
    keys = adict.keys()
    keys.sort()
    return [adict[key] for key in keys]

def all_files(root, patterns='*', single_level=True):
    """
    Iterate files under a given directory and for given patterns
    Arguments:
    o root - string, directory
    o patterns - string, patterns separated by ;
    o single_level - boolean, default: True
    """
    # Expand patterns from semicolon 
    patterns = patterns.split(';')
    for path, subdirs, files in os.walk(root):
        files.sort()
        for name in files:
            for pattern in patterns:
                if fnmatch.fnmatch(name, pattern):
                    yield os.path.join(path, name)
                    break
        if single_level:
            break

def get_dbpara(db_file):
    "Get database access information from a db parameter file"
    dbfile = check_file(db_file, PATH='/Users/yurong/src')
    dbfh = open(dbfile)
    for line in dbfh:
        line = line.rstrip()
        line_list = line.split('\t')
        if line_list[0] == 'dbhost':
            dbhost = line_list[1]
        elif line_list[0] == 'dbuser':
            dbuser = line_list[1]
        elif line_list[0] == 'dbpasswd':
            dbpasswd = line_list[1]
    dbfh.close()
    return (dbhost, dbuser, dbpasswd)

def search_file(filename, search_path):
    "Find file for given name and search path"
    candidate = os.path.join(search_path, filename)
    if os.path.isfile(candidate):
        return os.path.abspath(candidate)
    return None

def check_chr(chrN):
    "Check the chromosome name: chrN"
    if re.search('chr\w+', chrN):
        return chrN
    raise MethError('Wrong chromosome name %s' % chrN)

def check_file(filename, PATH=None):
    "Check if the file exits"
    if PATH is not None:
        candidate = os.path.join(PATH, filename)
    else:
        candidate = filename
    candidate = os.path.abspath(candidate)
    if os.path.isfile(candidate):
        return candidate
    raise MethError('Cannot find file %s' % candidate)

def count_file_lines(file):
    """
    Count lines of the file. 
    """
    count = -1
    for count, line in enumerate(open(file, 'rU')):
        pass
    count += 1
    return count
    
def finditer(text, pattern):
    "Find the pattern and return the iterator"
    pos = -1
    while True:
        pos = text.find(pattern, pos+1)
        if pos < 0:
            break
        yield pos

def contstr(attrbt, precision=None):
    """
    Convert other type to string.
    The main purpose of the function is to check None type.
    Arguments:
    o attrbt - attribute(s) of objects
    o precision - string of print format, default is '%.1f'
    type:
      o float
      o int (including LongType)
      o string
      o None
      o list
      o tuple
      o array
    Return values:
    o string of the attribute(s)
    """
    if type(attrbt) is NoneType or attrbt == []:
        return 'NULL'
    elif type(attrbt) is StringType:
        return attrbt
    elif type(attrbt) is IntType or type(attrbt) is LongType:
        return str(attrbt)
    elif type(attrbt) is FloatType:
        if precision == None:
            return "%.1f" % attrbt
        else:
            return precision % attrbt
    else: # list or other similar types
        str_attrbt = []
        for item in attrbt:
            if type(item) is FloatType:
                if precision == None:
                    str_attrbt.append("%.1f" % item)
                else:
                    str_attrbt.append(precision % item)
            elif type(item) is NoneType:
                str_attrbt.append('NULL')
            else:
                str_attrbt.append(str(item))
        return join(str_attrbt, "\t")

def read_seq(seqfh, format):
    "Get the sequence based on the format (fasta or genbank)"
    seq = ''
    if format == 'fasta':
        for line in seqfh:
            if re.search('^>', line):
                continue
            if re.search('\w', line):
                line = line.rstrip()
                seq += line
    else:
        raise MethError('Wrong sequence format')
    return seq.upper()


def get_sinfo(gfile):
    """
    Get table information for two groups
    Sample info file:
        # group 1
        table_name1   label1
        # group 2
        table_name2   label2
    """
    g1_tabinfo = {}  # sql table name: short sample name
    g2_tabinfo = {}
    gfh = open(gfile)
    for line in gfh:
        if re.search('^#', line) and re.search('group 1', line):
            tinfo = gfh.next().rstrip().split('\t')
            g1_tabinfo[tinfo[0]] = tinfo[1]
        elif re.search('^#', line) and re.search('group 2', line):
            tinfo = gfh.next().rstrip().split('\t')
            g2_tabinfo[tinfo[0]] = tinfo[1]
        else:
            continue
    gfh.close()
    return g1_tabinfo, g2_tabinfo

def get_sinfo_nogroup(gfile):
    """
    Samples are not divided into groups.
    Sample info file:
        table_name1   label1
        table_name2   label2
    Return:
    o tabinfo - dict
    """
    tabinfo = {}  # sql table name: short sample name
    gfh = open(gfile)
    for line in gfh:
        tinfo = line.rstrip().split('\t')
        tabinfo[tinfo[0]] = tinfo[1]
    gfh.close()
    return tabinfo

def cpgdensity_byseq(start, end, seq):
    """
    Compute CpG density for given coordinates and fasta file
    Input:
    o start - int, 0-based
    o end - int, 1-based
    o seq - Bio.Seq object
    Output:
    o CpG density - float
    """
    cpg_num = seq.count('CG', start, end)
    cpg_density = float(cpg_num * 2)/float(end - start)
    return cpg_density

