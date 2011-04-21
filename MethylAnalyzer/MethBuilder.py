# Classes used in all methylation landscape analyses

"""
MethBuilder collects classes and functions used in scripts for methylation analysis.
"""

import re
import sys
from math import floor
from numpy import mean
from pyfasta import Fasta

from MethylAnalyzer.MethError import MethError
from MethylAnalyzer.UtilityFuncs import finditer, EPSILON

# Global variables
MCRBC_SEQS = ['ACG', 'GCG', 'CGC', 'CGT']
# Corresponding enzymes in order: AciI, AciI_comp, HhaI, HpaII, HpyCH4V, BstUI
RE_SEQS = ['CCGC', 'GCGG', 'GCGC', 'CCGG', 'ACGT', 'CGCG']
# For CpG atlas table: column index vs. RE sequence
IDX2RE = {4:'CCGC', 5:'GCGG', 6:'GCGC', 7:'CCGG', 8:'ACGT', 9:'CGCG'}
# For CpG atlas table: column index vs. McrBC sequence
IDX2MCRBC = {10:'ACG', 11:'GCG', 12:'CGC', 13:'CGT'}

# Available Fasta sequences for hg18
HUMAN_CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
              'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
              'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM', 
              'chr1_random', 'chr2_random', 'chr3_random', 'chr4_random', 'chr5_random',
              'chr6_random', 'chr7_random', 'chr8_random', 'chr9_random', 'chr10_random',
              'chr11_random', 'chr13_random', 'chr15_random', 'chr16_random', 'chr17_random',
              'chr18_random', 'chr19_random', 'chr21_random', 'chr22_random', 'chrX_random',
              'chr5_h2_hap1', 'chr6_cox_hap1', 'chr6_qbl_hap2', 'chr22_h2_hap1']

HUMAN_MAIN_CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                   'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                   'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
              
# Available Fasta sequences for mm9
MOUSE_CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
              'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
              'chr18', 'chr19', 'chrX', 'chrY', 'chrM',
              'chr1_random', 'chr3_random', 'chr4_random', 'chr5_random', 'chr7_random',
              'chr8_random', 'chr9_random', 'chr13_random', 'chr16_random', 'chr17_random',
              'chrX_random', 'chrY_random', 'chrUn_random']

MOUSE_MAIN_CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                   'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                   'chr18', 'chr19', 'chrX', 'chrY']


# Classes
class Position:
    """
    Basic container class for nucleotides.
    This object has feature of coordinate.
    Arguments:
    o coordinate - int, coordinate for the nucleotide of interest
          e.g., the C in CpG, the first C in an RE recognition sequence CCGG
    """
    def __init__(self, coordinate):
        self.coordinate = coordinate

    def get_id(self):
        return self.coordinate

    def get_record(self, chr):
        """
        Provide genearl Position information
        Arguments:
        o chr - the specific chromosome that contains this position
        Record:
        chr, coordinate
        o coordinate - the first nucleotide in the recognition sequence (0-based)
        """
        record = '\t'.join([chr, str(self.coordinate)])
        return record


class PosDict:
    """
    The PosDict class defines a dict of Position objects. The Position objects here
    should come from the SAME chromosome. The analysis of methylation landscape is
    done by chromosome, so the coordinate of a Position object is a unique ID. The PosDict
    is indexed by coordinate or unique string id.
    Arguments:
    o label - string, the type of positions, e.g., 'CpG'
    o child_dict - dict, dict for child positions
    """
    def __init__(self, label):
        self.label = label
        self.child_dict = {}

    # Private methods
    def __len__(self):
        return len(self.child_dict)

    def __getitem__(self, id):
        "Return the child Position for given unique id"
        return self.child_dict[id]

    def __iter__(self):
        for value in self.child_dict.values():
            yield value        
            
    # Public methods
    def has_id(self, id):
        "Check if the PosDict has a Position with the given id"
        return self.child_dict.has_key(id)

    def add(self, position):
        "Add a child Position to the PosDict"
        pos_id = position.get_id()
        if pos_id in self.child_dict:
            raise MethError('%s position duplicated at %s' % (self.label, str(pos_id)))
        self.child_dict[pos_id] = position

    def delete(self, pos_id):
        "Delete a child position based on its id (coordinate)"
        del self.child_dict[pos_id]

    def sorted_iter(self):
        "Sort Position objects by coordinates"
        keys = self.child_dict.keys()
        keys.sort()
        for key in keys:
            yield self.child_dict[key]
            
    def sorted_coords(self):
        "Return sorted coordinates"
        keys = self.child_dict.keys()
        keys.sort()
        return keys


class CpGFull(Position):
    """
    The CpG class stores features of a CpG dinucleotide on a specific chromosome.
    Attributes:
    o coordinate - int, coordinate for the C in a CpG dinucleotide, 0-based
    o re_idx - int list, RE recognition sequence(s) that cover the CpG site
    o mcrbc_idx - int list, McrBC recognition sequence(s) that cover the CpG site
    o re_count - int, coverage for the RE library
    o mcrbc_count - int, coverage for the McrBC library
    o re_clean_count - int, coverage for the RE library after filtering
    o mcrbc_clean_count - int, coverage for the McrBC libary after filtering
    o reInterior - int, not include the two short ends of RE fragments [-> meth_value]
    o reEnd - int, short ends of RE fragments [-> unmeth_value]
    o mcrbcInterior - int, not include the two fuzzy ends of McrBC fragments [-> unmeth_value]
    """
    def __init__(self, coordinate):
        self.re_idx = []
        self.mcrbc_idx = []
        self.re_count = 0
        self.mcrbc_count = 0
        self.re_clean_count = 0
        self.mcrbc_clean_count = 0
        self.reInterior = 0
        self.reEnd = 0
        self.mcrbcInterior = 0
        Position.__init__(self, coordinate)

    # Public methods
    def dist2right(self, cpg):
        "Compute the distance to a CpG on the right"
        if cpg.coordinate <= self.coordinate:
            raise MethError('Wrong direction for two CpGs for distance calculation: CpG 1 <%d>, CpG 2 <%d>' % (self.coordinate, cpg.coordinate))
        self.dist = cpg.coordinate - self.coordinate

    def update_reseqs(self, idx):
        "Update RE recognition sequence using index to RE_SEQS"
        if idx not in self.re_idx:
            self.re_idx.append(idx)

    def update_mcrbcseqs(self, idx):
        "Update McrBC recognition sequence using index to MCRBC_SEQS"
        if idx not in self.mcrbc_idx:
            self.mcrbc_idx.append(idx)

    def is_re(self):
        "Check if the CpG overlaps with an RE site"
        if len(self.re_idx) > 0:
            return 1
        else:
            return 0

    def is_mcrbc(self):
        "Check if the CpG overlaps with an McrBC site"
        if len(self.mcrbc_idx) > 0:
            return 1
        else:
            return 0

    def get_recoverage(self):
        "The number of RE fragments hit the CpG"
        return self.re_count
    
    def get_mcrbccoverage(self):
        "The number of McrBC fragments hit the CpG"
        return self.mcrbc_count

    def get_reclean_coverage(self):
        return self.re_clean_count

    def get_mcrbcclean_coverage(self):
        return self.mcrbc_clean_count

    def get_coverage(self):
        "The total number of coverage"
        return self.re_count + self.mcrbc_count
    
    def get_methvalue(self):
        "The counts for methylation status"
        return self.reInterior

    def get_unmethvalue(self):
        "The counts for unmethylation status"
        return self.reEnd + self.mcrbcInterior

    def is_reInterior(self):
        "Check if the CpG is at the end of an RE fragment"
        if self.reInterior == 0:
            return 0
        else:
            return 1

    def is_reEnd(self):
        if self.reEnd == 0:
            return 0
        else:
            return 1

    def is_mcrbcInterior(self):
        if self.mcrbcInterior == 0:
            return 0
        else:
            return 1

    def cal_prob(self, p1_p2):
        """
        Estimate methylation probability for the CpG.
        Arguments:
        o p1_p2 - float, the ratio of sampling probability for RE and McrBC libraries
        Return:
        o p - float, methylation probability
        o n - float, the (theoretical) number of fragments for the CpG        
        """
        n1 = self.reInterior
        n2 = self.mcrbcInterior
        # Estimate p and n using the probability model
        self.n = n1 + n2
        try:
            self.p = n1/(n1 + p1_p2*n2)
        except ZeroDivisionError:
            raise MethError('ZeroDivisionError, Unable to estimate DNA methylation probability for CpG at position [ %d ]' % self.coordinate) 

    def get_record(self, chr):
        """
        Provide general CpG information
        Record:
        chr, coordinate, is_re, is_mcrbc, is_AciI, is_AciI_comp, is_HhaI, is_HpaII, is_HpyCH4V, is_BstUI,
        is_AGC, is_GCG, is_CGC, is_CGT
        return a string separated by \\t
        """
        is_re = str(self.is_re())
        is_mcrbc = str(self.is_mcrbc())
        re_states = ['0'] * 6
        for idx in self.re_idx:  # update positive hits for RE recognition sequence
            re_states[idx] = '1'
        re_states_str = '\t'.join(re_states)
        mcrbc_states = ['0'] * 4
        for idx in self.mcrbc_idx:  # update positive hits for McrBC recognition sequence
            mcrbc_states[idx] = '1'
        mcrbc_states_str = '\t'.join(mcrbc_states)
        record = '\t'.join([chr, str(self.coordinate), is_re, is_mcrbc, re_states_str, mcrbc_states_str])
        return record

    def get_methrecord(self, chr):
        """
        Methylation record for raw methylation data. The record can be used to rebuild CpG objects.
        Record:
           chr, coordinate, re_seqs, mcrbc_seqs, reInterior, reEnd, mcrbcInterior
           re_coverage, mcrbc_coverage, re_clean_coverage, mcrbc_clean_coverage,
           is_re, is_mcrbc
        """
        re_seqs = mcrbc_seqs = 'NULL'
        if self.re_idx:
            re_seqs = ','.join([RE_SEQS[idx] for idx in self.re_idx])
        if self.mcrbc_idx:
            mcrbc_seqs = ','.join([MCRBC_SEQS[idx] for idx in self.mcrbc_idx])
        reInterior = str(self.reInterior)
        reEnd = str(self.reEnd)
        mcrbcInterior = str(self.mcrbcInterior)
        re_coverage = str(self.re_count)
        re_clean_coverage = str(self.re_clean_count)
        mcrbc_coverage = str(self.mcrbc_count)
        mcrbc_clean_coverage = str(self.mcrbc_clean_count)
        is_re = str(self.is_re())
        is_mcrbc = str(self.is_mcrbc())
        record = '\t'.join([chr, str(self.coordinate), re_seqs, mcrbc_seqs, reInterior,
                            reEnd, mcrbcInterior, re_coverage, mcrbc_coverage,
                            re_clean_coverage, mcrbc_clean_coverage, is_re, is_mcrbc])
        return record

    def get_tablerecord(self, chr):
        """
        Generate the record for a database table.
        Record:
           chr, coordinate, is_remix, reInterior, reEnd, mcrbcInterior, re_clean_coverage,
           mcrbc_clean_coverage
        o is_remix - the CpG is in an RE site and it is covered by both reEnd and reInterior
        """
        if self.is_reEnd() and self.is_reInterior() and self.is_re():
            is_remix = '1'
        else:
            is_remix = '0'
        record = map(str, [chr, self.coordinate, is_remix, self.reInterior, self.reEnd,
                           self.mcrbcInterior, self.re_clean_count, self.mcrbc_clean_count])
        record = '\t'.join(record)
        return record
        
    def get_coverrecord(self, chr):
        """
        Generate the record for RE/McrBC coverage. A CpG is covered by an RE/McrBC fragment if it
        resides within the fragment.
        Record:
        chr, coordinate, RE raw, RE filtered, McrBC raw, McrBC filtered
        """
        record = '\t'.join(map(str, [chr, self.coordinate, self.re_count, self.re_clean_count,
                                     self.mcrbc_count, self.mcrbc_clean_count]))
        return record


class Fragment:
    """
    The Fragment class is used to process fragments parsed from *.mates files (see FragmentParser).

    Attributes:
    o id - string, unique id of the fragment, original id + '-' + left_coord + '-' right_coord
    o left_coord - int, left coordinate (0-based), start coordinate, smaller number
    o right_coord - int, right coordinate (1-based), end coordinate, larger number
    o strand - char, + or -
    o readlib - string, RE or McrBC
    o filter - int, status after filtering
           1: pass the filtering with RE/McrBC site at both ends;
           2: pass the filtering with RE/McrBC site at only one end;
           3: fail, no RE/McrBC site at the two ends;
           4: fail, no RE/McrBC site in the interior region of the fragment
           
    Attributes genearted by function
    o leftsub_in - int, the coordinate inside of the left boundary
    o leftsub_out - int, the coordinate outside of the left boundary
    o rightsub_in - int, the coordinate inside of the right boundary
    o rightsub_out - int, the coordinate outside of the right boundary
    """
    def __init__(self, id, left_coord, right_coord, strand, readlib):
        self.id = id
        self.left_coord = left_coord
        self.right_coord = right_coord
        self.strand = strand
        self.readlib = readlib
        self.filter = None

    # Private methods
    def __len__(self):
        "Return the length of the fragment"
        return self.right_coord - self.left_coord

    # Public methods
    def get_id(self):
        "id = mates_id + '-' + left_coord + '-' + right_coord"
        return self.id

    def get_matesid(self):
        "Get original id provided by the mates file"
        return self.id.split('-')[0]

    def get_leftout(self, outlen):
        """
        Define the outside coordinate of the fuzzy end for the start boundary.
        o outlen - the parameter (sequence length) to generate the fuzzy end
        """
        return self.left_coord - outlen

    def get_leftin(self, inlen, off):
        """
        Define the inside coordinate of the fuzzy end for the start boundary.
        o inlen - int, the parameter (sequence length) to generate the fuzzy end
        o off - int, based on the enzyme recognition sequence, minus this number to
            search the complete enzyme recognition seq within the fuzzy end
        """
        return self.left_coord + inlen - off

    def get_rightin(self, inlen):
        """
        Define the inside coordinate of the fuzzy end for the end boundary.
        Since all the start coordinates are 0-based, i.e., coordinates of Position objects
        (e.g., CpG, RE site, McrBC site), the end coordinate should be -1 in order
        to reset it to 0-based.
        o inlen - int, the parameter (sequence length) to generate the fuzzy end
        """
        return self.right_coord - 1 - inlen

    def get_rightout(self, outlen, off):
        """
        Define the outside coordinate of the subsequence for the end boundary.
        Since all the start coordinates are 0-based, i.e., coordinates of Position objects
        (e.g., CpG, RE site, McrBC site), the end coordinate should be -1 in order to
        reset it to 0-based.
        o outlen - int, the parameter (sequence length) to generate the subsequence
        o off - int, based on the enzyme recognition sequence, minus this number to
                search the complete enzyme recognition seq within the subsequence
        """
        return self.right_coord - 1 + outlen - off

    def get_bedrecord(self, chr):
        """
        Provide general fragment information including in 'bed' format.
        The records can be viewed in UCSC Genome Browser.
        Record:
           chr, start coordinate, end coordinate, name, score (use 1 here), strand,
           thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
        """
        score = '1'
        itemRgb = '0'
        blockCount = '2'
        blockSizes = '25,25,'
        blockStarts = '0,' + str(self.right_coord - self.left_coord - 25)
        mates_id = self.get_matesid()  # original id
        record = '\t'.join([chr, str(self.left_coord), str(self.right_coord), mates_id, score, 
                            self.strand, str(self.left_coord), str(self.right_coord), itemRgb,
                            blockCount, blockSizes, blockStarts])
        return record

    def get_filterrecord(self, chr):
        """
        Provide detailed information after RE/McrBC site filtering.
        Record:
           chr, start coordinate, end coordinate, name, strand, readlib, filter
        """
        mates_id = self.get_matesid()  # original id
        record = '\t'.join([chr, str(self.left_coord), str(self.right_coord), mates_id,
                            self.strand, self.readlib, str(self.filter)])
        return record


class SiteParser:
    """
    Parse all CpG-related sites that are needed in the methylation profiling,
    including CpG and RE/McrBC recognition sites.

    Attributes:
    o chr - string, chromosome
    o sites - PosDict, sites for CpG or RE or McrBC
    
    Additional attributes:
    o seqfile - string, chromosome fasta sequence filename
    o seq - string, the sequence
    o sitefile - string, filename for parsed sites
    """
    def __init__(self, chr):
        self.chr = chr
        self.sites = PosDict('Site')

    # Private methods
    def _parse_cpgs(self, sitefh):
        for line in sitefh:
            line_list = line.rstrip().split('\t')
            # Get values
            coord = int(line_list[1])
            re_idx = []
            mcrbc_idx = []
            for i in range(4, 10):  # for six RE recognition sequences
                if int(line_list[i]) == 1:
                    re_idx.append(i-4)  # convert to index of RE_SEQS
            for i in range(10, 14):  # for four McrBC recognition sequences
                if int(line_list[i]) == 1:
                    mcrbc_idx.append(i-10)  # convert to index of MCRBC_SEQS
            # Initiate CpG
            cpg = CpGFull(coord)
            cpg.re_idx = re_idx
            cpg.mcrbc_idx = mcrbc_idx
            try:
                self.sites.add(cpg)
            except MethError, e:
                print >> sys.stderr, 'MethError: ', e.value
                sys.exit(2)
                
    def _parse_enzymes(self, sitefh):
        for line in sitefh:
            line_list = line.rstrip().split('\t')
            coord = int(line_list[1])
            enzyme = Position(coord)
            try:
                self.sites.add(enzyme)
            except MethError, e:
                print >> sys.stderr, 'MethError: ', e.value
                sys.exit(2)
    
    # Public methods
    def get_sites(self):
        return self.sites

    def parse_seq(self, seqfile):
        """
        Scan the given sequence (usually one chromosome) to generate CpG/RE/McrBC sites
        Note: fasta file has head(s) same as chromosome name(s)
        """
        fa_record = Fasta(seqfile)
        seq = str(fa_record[self.chr]).upper()
        # 1. Parse CpG sites, use the coordinate of C
        for coord in finditer(seq, 'CG'):
            cpg = CpGFull(coord)
            try:
                self.sites.add(cpg)
            except MethError, e:
                print >> sys.stderr, 'MethError: ', e.value
                sys.exit(2)
        # 2. Parse RE sites, use the coordinate of the first nucleotide
        #    in the recognition sequence
        self.re_sites = PosDict('RE')
        for re_seq in RE_SEQS:
            idx = RE_SEQS.index(re_seq)
            for coord in finditer(seq, re_seq):
                re = Position(coord)
                try:
                    self.re_sites.add(re)
                except MethError, e:
                    print e.value
                    sys.exit(2)
                if re_seq == 'CGCG':
                    # Update the CpG object in the RE site at coord
                    self.sites[coord].update_reseqs(idx)
                    # Update the CpG object in the RE site at coord + 2
                    self.sites[coord + 2].update_reseqs(idx)                
                else:
                    # Update the CpG object in the RE site at coord + 1
                    self.sites[coord + 1].update_reseqs(idx)
        # 3. Parse McrBC sites, use the coordinate of the first nucleotide
        #    in the recognition sequence
        self.mcrbc_sites = PosDict('McrBC')
        for mcrbc_seq in MCRBC_SEQS:
            idx = MCRBC_SEQS.index(mcrbc_seq)
            for coord in finditer(seq, mcrbc_seq):
                mcrbc = Position(coord)
                try:
                    self.mcrbc_sites.add(mcrbc)
                except MethError, e:
                    print e.value
                    sys.exit(2)
                if mcrbc_seq == 'ACG' or mcrbc_seq == 'GCG':  # for 'RCG'
                    # Update the CpG object in the McrBC site at coord + 1
                    self.sites[coord + 1].update_mcrbcseqs(idx)
                else:  # for 'CGY'
                    # Update the CpG object in the McrBC site at coord
                    self.sites[coord].update_mcrbcseqs(idx)

    def parse_sites(self, sitefile, label):
        """
        Build sites from the sitefile (stored in the annotation directory)
        Input
        o sitefile - string, filename
        o label - string, CpG/RE/McrBC
        
        CpG annotation file format:
           chr, coordinate, is_re, is_mcrbc, is_AciI, is_AciI_comp, is_HhaI
           is_HpaII, is_HpyCH4V, is_BstUI, is_AGC, is_GCG, is_CGC, is_CGT
        RE/McrBC annotation format:
           chr  coordinate
        """
        sitefh = open(sitefile)
        # Parse CpGs
        if label == 'CpG':
            self._parse_cpgs(sitefh)
        # Parse RE/McrBC sites
        elif label == 'RE' or label == 'McrBC':
            self._parse_enzymes(sitefh)
        sitefh.close()

    def cpgdensity_bycpg(self, winlen):
        """
        Compute CpG density by CpG-centered windows. The method should be applied after
        'parse_sites', i.e., this method needs CpG site information to compute CpG density.
        Input:
        o winlen - int, window length centered by the current CpG
        """
        if len(self.cpg_sites) == 0:
            raise MethError("Need to parse CpG sites first")
        for cpg in self.cpg_sites.sorted_iter():
            cpg_num = 0
            position = cpg.coordinate
            left_bound = position - winlen/2
            right_bound = position + winlen/2
            # Search CpGs within the defined window centered by the current CpG
            for pos in range(left_bound, right_bound):
                if self.cpg_sites.has_id(pos):
                    cpg_num += 1
            cpg_density = float(2 * cpg_num)/float(right_bound - left_bound)
            cpg.density = cpg_density
            
    def oeratio_bycpg(self, winlen, fafile):
        """
        Compute CpG O/E ratio by CpG-centered windows. The method should be applied after
        'parse_sites', i.e., this method needs CpG site information to compute CpG O/E ratio.
        Input:
        o winlen - int, window length centered by the current CpG
        o fafile - file, Fasta file
        """
        # Parse Fasta file
        fa_record = Fasta(fafile)
        seq = str(fa_record[self.chr]).upper()
        seqlen = len(seq)
        # Compute O/E ratio
        for cpg in self.cpg_sites.sorted_iter():
            position = cpg.coordinate
            start = position - winlen/2  # 0-based
            end = position + winlen/2  # 1-based
            if start < 0:
                start = 0
                end = winlen
            if end > seqlen:
                end = seqlen
                start = seqlen - winlen
            cpg_num = seq.count('CG', start, end)
            c_num = seq.count('C', start, end)
            g_num = seq.count('G', start, end)
            oe_ratio = float(cpg_num * winlen)/float(c_num * g_num)
            cpg.oe_ratio = oe_ratio

    def cpgdensity_byregion(self, start, end):
        """
        Compute CpG density by CpG-centered windows. The method should be applied after
        'parse_sites', i.e., this method needs CpG site information to compute CpG density.
        Input:
        o start - int, 0-based coordinate
        o end - int, 1-based coordinate
        Output
        o cpg_density - float
        """
        if len(self.cpg_sites) == 0:
            raise MethError("Need to parse CpG sites first")
        # Search CpGs within the defined window centered by the current CpG
        for pos in range(start, end):
            if self.cpg_sites.has_id(pos):
                cpg_num += 1
        cpg_density = float(2*cpg_num)/float(end-start)
        return cpg_density
    

class FragmentParser:
    """
    Parse paired-end reads to obtain coordinates of methylaed/unmethylated fragments.
    Attributes:
    o matesfile - string, paired-end file name
    o chr_mapping - dict, map chromosome indeces to chromosome names, such as {num: chr},
                    employed by SOLiD aglinment 
    o fragfile - string, file name for the parsed fragments (output)
    o readlib - string, RE or McrBC
    """
    def __init__(self, CHRMAP=None):
        if CHRMAP is not None:
            self.chr_mapping = self._get_mapping(CHRMAP)

    # Private methods
    def _get_mapping(self, cmapfile):
        """
        Get the chromosome mapping.
        Format of the cmap file:
           index(int)   chromosome   dir   ...
        """
        chr_mapping = {}
        mapfh = open(cmapfile)
        for line in mapfh:
            if re.search('^\d+\t', line):
                line_list = line.rstrip().split('\t')
                chr_mapping[line_list[0]] = 'chr' + line_list[1]
        mapfh.close()
        return chr_mapping

    # Public methods
    def parse_mates(self, matesfile, readlib, ERRORFILE='frag_errors'):
        """
        Read the mates file and yield fragment iteratively.
        Check the document of SOLiD data format for further information of *.mates files.
        Note: coordinates follow the UCSC Genome Browser tridition,
              0-based for start and 1-based for end.
        """
        matesfh = open(matesfile)
        errorfh = open(ERRORFILE, 'a')
        for line in matesfh:
            # Skip the header lines
            if re.search('^#', line):
                continue
            line_list = line.rstrip().split('\t')
            # Parse only fragments in the category of 'AAA'
            if line_list[-1] != 'AAA':
                continue
            ftag_len = len(line_list[1]) - 1
            chr = self.chr_mapping[line_list[6]]
            if int(line_list[8]) < 0:
                strand = '-'
                left_coord = abs(int(line_list[8])) - ftag_len + 1  # 0-based
                right_coord = abs(int(line_list[9])) + 1  # 1-based
            else:
                strand = '+'
                left_coord = int(line_list[9])  # 0-based
                right_coord = int(line_list[8]) + ftag_len  # 1-based
            # Check fragments with errors
            if left_coord >= right_coord:
                errorfh.write('Fragment aligned wrong [%s]: %s, %s strand, %d-%d\n' % \
                              (line_list[0], chr, strand, left_coord, right_coord))
                continue
            # The original id provided by mates file may not be unique
            # Redefine: fragment id = original id + '-' + start + '-' + end
            id = line_list[0] + '-' + str(left_coord) + '-' + str(right_coord)
            fragment = Fragment(id, left_coord, right_coord, strand, readlib)
            yield chr, fragment
        matesfh.close()
        errorfh.close()

    def parse_bedfrags(self, fraglines, readlib):
        """
        Generate fragments (BED format) generated by parse_mates.
        By doing that, the program can split large file into several subsets
        in order to avoid memory overflow.
        BED format (Note: there are headers for the bed format): 
           chr, start coordinate, end coordinate, name, score (use 1 here), strand,
           thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts           
        """
        for line in fraglines:
            # Skip header lines
            if not re.search('^chr\w+', line):
                continue
            line_list = line.rstrip().split('\t')
            left_coord = int(line_list[1])
            right_coord = int(line_list[2])
            id = line_list[3]
            strand = line_list[5]
            fragment = Fragment(id, left_coord, right_coord, strand, readlib)
            yield fragment


class REFilter:
    """
    Filter RE fragments based on a set of criteria.

    Attributes:
    o re_frags - list of RE fragments
    o cpg_sites - PosDict, CpG sites
    o re_sites - PosDict, RE sites
    o outlen - int, the length outside of the 5'- and 3'- fragment ends to serve as 'fuzzy' ends
    o inlen - int, the length inside of the 5'- and 3'- fragment ends to serve as 'fuzzy' ends

    Procedure:
    1. Search RE sites in the 5'- and 3'- fuzzy ends defined by outlen and inlen
       Note: RE sites should reside within the fuzzy ends, don't consider partially
       overlapping sites
    1a. Get the left fuzzy end (i.e., 5'-end referring to the + strand)
        out: outside of the boundary, in: inside of the boundary
        Use 0-based convention for both coordinates
    1b. Get the right fuzzy end, same as 1a
    1c. Get the sequence in the middle without the fuzzy ends
    2. Check if at least one fuzzy end has a valid RE recognition site
    3. Update fragment status and CpG status
    4. Update the number of fragment coverage for the CpG
       include 'total' converage and 'clean' coverage after filtering    
    """
    def __init__(self, re_frags, cpg_sites, re_sites, outlen, inlen):
        self.re_frags = re_frags
        self.cpg_sites = cpg_sites
        self.re_sites = re_sites
        self.outlen = outlen
        self.inlen = inlen

    # Private methods
    def _search_site(self, left, right):
        "Search RE sites in the fuzzy end."
        for coord in range(left, right):
            if self.re_sites.has_id(coord):
                return True
        return False
                
    def _update_cpgs(self, left_out, left_in, right_in, right_out, left_re, right_re, middle_re):
        """
        Update methylation status for CpGs in the RE fragment.
        The CpGs in the fuzzy ends and overlap with RE sites are assigned to reEnd.
        The CpGs in the middle are assigned to reInterior.
        Note: for CpGs at the fuzzy ends, only update the CpGs that overlap with RE sites
        """
        # For the left fuzzy end
        # Update CpGs that overlap with RE sites
        if left_re:
            for coord in range(left_out, left_in):
                if self.cpg_sites.has_id(coord) and self.cpg_sites[coord].is_re():
                    self.cpg_sites[coord].reEnd += 1
        # For the right fuzzy end
        # Update CpGs that overlap with RE sites
        if right_re:
            for coord in range(right_in, right_out):
                if self.cpg_sites.has_id(coord) and self.cpg_sites[coord].is_re():
                    self.cpg_sites[coord].reEnd += 1
        # For the interior subsequence [left_in, right_in]
        # Update only if there is at least one RE site in the interior subsequence
        if middle_re:
            for coord in range(left_in, right_in):
                if self.cpg_sites.has_id(coord):
                    self.cpg_sites[coord].reInterior += 1
        
    def _update_fragcount(self, fragment):
        """
        Search all CpGs within the fragment.
        If the C is on the right boundary (i.e., the last nucleotide is C), it would not be updated.
        """
        for coord in range(fragment.left_coord, fragment.right_coord-1):  # make sure only intact CpGs would be found in a fragment
            if self.cpg_sites.has_id(coord):
                self.cpg_sites[coord].re_count += 1
                if fragment.filter == 1 or fragment.filter == 2:
                    self.cpg_sites[coord].re_clean_count += 1

    # Public methods
    def scan(self):
        "Filter RE fragments."
        for fragment in self.re_frags:
            # 1. Search RE sites in the two fuzzy ends
            #    Note: RE sites should reside within the subseq, don't consider partially overlapping sites
            # 1a. Left fuzzy end, out: outside of the boundary, in: inside of the boundary
            #    Both coordinates use 0-based convention
            left_out = fragment.get_leftout(self.outlen)
            left_in = fragment.get_leftin(self.inlen, 2)  # off=2, the last possible RE site: NNN(S) S is the boundary
            left_re = self._search_site(left_out, left_in)
            # 1b. Right fuzzy end
            right_in = fragment.get_rightin(self.inlen)
            right_out = fragment.get_rightout(self.outlen, 2)  # off = 2
            right_re = self._search_site(right_in, right_out)
            # 1c. Sequence in the middle without the fuzzy ends
            middle_re = self._search_site(left_in, right_in)
            # 2. Update fragment status and CpG status
            if left_re or right_re:
                self._update_cpgs(left_out, left_in, right_in, right_out, left_re, right_re, middle_re)
                if middle_re:
                    if left_re and right_re:
                        fragment.filter = 1  # passed: both ends have RE sites
                    else:
                        fragment.filter = 2  # passed: only one end has RE sites
                else:
                    fragment.filter = 4  # failed: no RE sites in the middle
            else:
                fragment.filter = 3  # failed: no RE sites at the two ends
            # 3. Update the number of fragment coverage for the CpG
            # include 'total' converage and 'clean' coverage after filtering
            self._update_fragcount(fragment)

    def get_passed(self):
        "Get fragments that pass the filtering step."
        frag_passed = []
        for fragment in self.re_frags:
            if fragment.filter == 1 or fragment.filter == 2:
                frag_passed.append(fragment)
        return frag_passed

    def get_failed(self):
        "Get the failed fragments in two categories: end and middle."
        frag_failed_end = []
        frag_failed_mid = []
        for fragment in self.re_frags:
            if fragment.filter == 3:
                frag_failed_end.append(fragment)
            elif fragment.filter == 4:
                frag_failed_mid.append(fragment)
        return frag_failed_end, frag_failed_mid

    def get_cpgs(self):
        return self.cpg_sites


class MCRBCFilter:
    """
    Filter McrBC fragments based on a set of criteria.

    Attributes:
    o mcrbc_frags - PosDict, McrBC fragments
    o cpg_sites - PosDict, CpG sites
    o mcrbc_sites - PosDict, McrBC sites
    o winlen - int, the length inside or outside of the 5'- and 3'-fragment ends to
               search for McrBC half sites, default 500bp
    o cut_dist - int, the distance from the end, because the cleavage position occurs about
                 30bp from the McrBC half site (REF: Panne et al, JMB, 290:49-60, 1999)
                 default value: 50bp
    o min_dist - int, the minimum distance between two McrBC half sites, default 40bp
    o max_dist - int, the maximum distance between two McrBC half sites, default 500bp
                 (NOTE: generally the distance between two half sites is from 40bp to 500bp)
    o fuzzylen - int, the length used to determine the interior segment of an McrBC fragment, default 50bp

    Procedure:
    1. Search the two McrBC half sites in the two fragment ends. A search region is split into two
       segments by the actual fragment end. Each segment should reside one McrBC half site.
    1a. For the left end: (1) find at least one McrBC half site in each segment;
        (2) the two half sites should satisfy the minimum and maximum distance criteria;
        (3) and one of the half sites should reside within the cut distance from the fragment end.
    1b. For the right end, run same search
    1c. Search a valid McrBC site (two half sites) in the middle, i.e., minus the fuzzy ends
        the fuzzy ends are defined by fuzzylen
    2. Update fragment status and CpG status
    3. Update the number of fragment coverage for the CpG
       include 'total' converage and 'clean' coverage after filtering
    """
    def __init__(self, mcrbc_frags, cpg_sites, mcrbc_sites, winlen, cut_dist, min_dist, max_dist, fuzzylen):
        self.mcrbc_frags = mcrbc_frags
        self.cpg_sites = cpg_sites
        self.mcrbc_sites = mcrbc_sites
        self.winlen = winlen
        self.cut_dist = cut_dist
        self.min_dist = min_dist
        self.max_dist = max_dist
        self.fuzzylen = fuzzylen

    # Private methods
    def _search_midsite(self, left, right):
        """
        Search the McrBC site in the middle subsequence.
        Return True if there is at least one pair of half sites and two sites > the mininum distance
        and < the maximum distance
        """
        sites = []
        # Scan from left
        for coord in range(left, right):
            if self.mcrbc_sites.has_id(coord):
                sites.append(coord)
        # Check the distance criteria
        for m in range(0, len(sites)-1):
            for n in range(m+1, len(sites)):
                if abs(sites[n] - sites[m]) > self.min_dist and abs(sites[n] - sites[m]) < self.max_dist:
                    return True
        return False
        
    def _search_endsite(self, left, right, boundary):
        """
        Search the McrBC half site in [left, boundary] and [boundary, right].
        Procedure:
        1. find all McrBC half sites in two segments [left, boundary] and [boundary, right]
        2. get the near and far half sites in terms of the boundary
        3. check if the near site is within the cut distance from the boundary
        4. check if the distance between the satisfied near site at one side and
           the far site at the other side is > the minimum distance
        """
        # 1. Find all McrBC half sites
        left_sites = []
        for coord in range(left, boundary):
            if self.mcrbc_sites.has_id(coord):
                left_sites.append(coord)
        right_sites = []
        for coord in range(boundary, right):
            if self.mcrbc_sites.has_id(coord):
                right_sites.append(coord)
        # 2. Get near and far sites
        # near: the closest site to the boundary
        # far: the farthest site to the boundary
        if left_sites == [] or right_sites == []:
            return False
        left_near = max(left_sites)
        left_far = min(left_sites)
        right_near = min(right_sites)
        right_far = max(right_sites)
        # 3. Check if at least one McrBC half site is within the cut_dist
        #    and check the minimum distance between this near site
        #    and the far site on the other side
        if left_near >= boundary - self.cut_dist and abs(left_near - right_far) > self.min_dist:
            return True
        if right_near <= boundary + self.cut_dist - 2 and abs(right_near - left_far) > self.min_dist:
            return True
        return False

    def _update_cpgs(self, left_in, right_in, middle_mcrbc):
        """
        Update methylation status for CpGs in the McrBC fragment.
        The CpGs in the inside subsequences are assigned to mcrbcInterior.
        """
        # For the interior subsequence [left_in, right_in]
        # Update only if there is at least an McrBC site in the interior subsequence
        if middle_mcrbc:
            for coord in range(left_in, right_in):
                if self.cpg_sites.has_id(coord):
                    self.cpg_sites[coord].mcrbcInterior += 1
        
    def _update_fragcount(self, fragment):
        """
        Update all CpGs within the fragment.
        If the C is on the right boundary (i.e., the last nt is C), it would not be updated.
        """
        for coord in range(fragment.left_coord, fragment.right_coord-1):  # make sure only intact CpGs would be found in a fragment
            if self.cpg_sites.has_id(coord):
                self.cpg_sites[coord].mcrbc_count += 1
                if fragment.filter == 1 or fragment.filter == 2:
                    self.cpg_sites[coord].mcrbc_clean_count += 1

    # Public methods
    def scan(self):
        "Filter McrBC fragments."
        for fragment in self.mcrbc_frags:
            # 1. Search the two McrBC half sites in the two fragment ends
            # 1a. For the left boundary
            leftsearch_out = fragment.get_leftout(self.winlen)
            leftsearch_in = fragment.get_leftin(self.winlen, 1)  # off = 1
            left_mcrbc = self._search_endsite(leftsearch_out, leftsearch_in, fragment.left_coord)
            # 1b. For the right boundary
            rightsearch_in = fragment.get_rightin(self.winlen)
            rightsearch_out = fragment.get_rightout(self.winlen, 1)  # off = 1
            right_mcrbc = self._search_endsite(rightsearch_in, rightsearch_out, fragment.right_coord)
            # 1c. Sequence in the middle without the fuzzy ends
            # the fuzzy ends are defined by fuzzylen
            left_out = fragment.get_leftout(self.fuzzylen)  #### to del
            left_in = fragment.get_leftin(self.fuzzylen, 1)  # off = 1
            right_in = fragment.get_rightin(self.fuzzylen)
            right_out = fragment.get_rightout(self.fuzzylen, 1)  # off = 1, to del
            middle_mcrbc = self._search_midsite(left_in, right_in)
            # 2. Update fragment status and CpG status
            if left_mcrbc or right_mcrbc:
                self._update_cpgs(left_in, right_in, middle_mcrbc)
                if middle_mcrbc:
                    if left_mcrbc and right_mcrbc:
                        fragment.filter = 1  # passed: both ends have McrBC sites
                    else:
                        fragment.filter = 2  # passed: only one end has McrBC sites
                else:
                    fragment.filter = 4  # failed: no McrBC sites in the middle
            else:
                fragment.filter = 3  # failed: no McrBC sites in the fuzzy ends 
            # 3. Update the number of fragment coverage for the CpG
            # include 'total' converage and 'clean' coverage after filtering
            self._update_fragcount(fragment)
                
    def get_passed(self):
        "Get fragments that pass the filtering step."
        frag_passed = []
        for fragment in self.mcrbc_frags:
            if fragment.filter == 1 or fragment.filter == 2:
                frag_passed.append(fragment)
        return frag_passed

    def get_failed(self):
        "Get the failed fragments in two categories: end and middle"
        frag_failed_end = []
        frag_failed_mid = []
        for fragment in self.mcrbc_frags:
            if fragment.filter == 3:
                frag_failed_end.append(fragment)
            elif fragment.filter == 4:
                frag_failed_mid.append(fragment)
        return frag_failed_end, frag_failed_mid

    def get_cpgs(self):
        return self.cpg_sites


class MethStateParser:
    """
    Determine methylation status of CpG sites.
    Attributes:
    o chr - string, chromosome, e.g., chr11
    o chr_len - int, length of the chromosome
    o space - int, space used to split the chromosome
    o cpg_sites - PosDict for CpGs
    CpGs with methvalue + unmethvalue > 0
                  CpGs have been processed for the RE/McrBC filtering
                  use 'set_cpgsites' or 'rebuild_cpgs' to assign the cpg_sites
    """
    def __init__(self, chr_len, space):
        self.chr_len = chr_len
        self.space = space

    # Private methods
    def _cal_p_bar(self):
        """
        Calculate p_bar based on newly estimated methylation probabilities.
        p_bar = mean([pi, ...])
        pi = mean(p of CpGs in segment i)
        """
        pis = []
        for start in range(0, self.chr_len, self.space):
            end = start + self.space  # 1-based
            inseg_ps = []  # initiate the list for each segment
            for pos in range(start, end-1):  # make sure only intact CpGs will be found in a segment
                if self.cpg_sites.has_id(pos):
                    inseg_ps.append(self.cpg_sites[pos].p)
            if len(inseg_ps) > 0:  # the segment contains at least one CpG
                pave = mean(inseg_ps)
                pis.append(pave)
        return mean(pis)
    
    def _parse_methline(self, line):
        "Parse one line from the methylation data file"
        line_list = line.rstrip().split('\t')
        coord = int(line_list[1])
        # Methylation status counts
        reInterior = int(line_list[4])
        reEnd = int(line_list[5])
        mcrbcInterior = int(line_list[6])
        # Indexes to recognition sequences
        if line_list[2] != 'NULL':
            re_idx = [RE_SEQS.index(seq) for seq in line_list[2].split(',')]
        else:
            re_idx = []                
        if line_list[3] != 'NULL':
            mcrbc_idx = [MCRBC_SEQS.index(seq) for seq in line_list[3].split(',')]
        else:
            mcrbc_idx = []
        # Fragment counts
        re_count = int(line_list[7])
        mcrbc_count = int(line_list[8])
        re_clean_count = int(line_list[9])
        mcrbc_clean_count = int(line_list[10])
        return coord, re_idx, mcrbc_idx, re_count, mcrbc_count, re_clean_count, mcrbc_clean_count, reInterior, reEnd, mcrbcInterior

    def _sumn_in_segments(self):
        "Sum of average RE and McrBC read coverage of genomic segments defined by the parser."
        n1_list = []
        n2_list = []
        for start in range(0, self.chr_len, self.space):
            end = start + self.space  # 1-based
            inseg_n1 = []  # initiate the list for each segment
            inseg_n2 = []
            for pos in range(start, end-1):  # make sure only intact CpGs will be found in a segment
                if self.cpg_sites.has_id(pos):
                    inseg_n1.append(self.cpg_sites[pos].reInterior)
                    inseg_n2.append(self.cpg_sites[pos].mcrbcInterior)
            if len(inseg_n1) > 0:  # the segment contains at least one CpG
                n1_list.append(mean(inseg_n1))
                n2_list.append(mean(inseg_n2))
        return sum(n1_list), sum(n2_list)
    
    def _get_sampling_paras(self, p_bar, sum_n1, sum_n2):
        """
        Estimate sampling probabilities for each chromosome.
        A chromosome is split into equally spaced segments with length of 1000bp.
        The total number of segments is N which is independent of the RE or McrBC library.
        For a segment i,
           n1,i = mean(n1's in segment i)   n1 = reInterior     
           n2,i = mean(n2's in segment i)   n2 = mcrbcInterior  
        a1 = sum(n1,i)/p_bar
        a2 = sum(n2,i)/(1-p_bar)
        p1/p2 (lamda) ~ a1/a2

        Note: p_bar is determined by LUMA assay.
        """
        a1 = sum_n1/p_bar
        a2 = sum_n2/(1-p_bar)
        try:
            p1_p2 = a1/a2
        except ZeroDivisionError:
            raise MethError('ZeroDivisionError, Unable to estimate sampling probability')
        # Check lamda value for the extremely biased sampling between RE and McrBC
        if p1_p2 > 10.0 or p1_p2 < 0.1:
            raise MethError('Extremely biased sampling between RE and McrBC, lamda = %.2f' % p1_p2)
        return p1_p2
        

    # Public methods
    def set_cpgsites(self, cpg_sites):
        "CpG sites parsed from the CpG atlas file"
        self.cpg_sites = cpg_sites

    def get_cpgs(self):
        return self.cpg_sites
        
    def parse_methdata(self, methfile):
        """
        Update count information to the pre-built CpG sites.
        Note that CpG sites should have been initialized by the method
        'set_cpgsites'
        """
        methfh = open(methfile)
        for line in methfh:
            coord, re_idx, mcrbc_idx, re_count, mcrbc_count, re_clean_count, mcrbc_clean_count, reInterior, reEnd, mcrbcInterior = self._parse_methline(line)
            # Skip the line if the methylation value is equal to 0
            if reInterior + mcrbcInterior == 0:
                continue
            cpg = self.cpg_sites[coord]
            # Update RE/McrBC counts
            cpg.re_count += re_count
            cpg.mcrbc_count += mcrbc_count
            cpg.re_clean_count += re_clean_count
            cpg.mcrbc_clean_count += mcrbc_clean_count
            cpg.reInterior += reInterior
            cpg.reEnd += reEnd
            cpg.mcrbcInterior += mcrbcInterior
        methfh.close()

    def rebuild_cpgs(self, methfile):
        """
        Rebuild CpG objects from the methfile generated by filter.py.
        Methylation records:
           chr, coordinate, re_seqs, mcrbc_seqs, reInterior, reEnd, mcrbcInterior
           re_coverage, mcrbc_coverage, re_clean_coverage, mcrbc_clean_coverage,
           is_re, is_mcrbc
        """
        self.cpg_sites = PosDict('CpG')
        cpgfh = open(methfile)
        for line in cpgfh:
            coord, re_idx, mcrbc_idx, re_count, mcrbc_count, re_clean_count, mcrbc_clean_count, reInterior, reEnd, mcrbcInterior = self._parse_methline(line)
            # Skip "uncovered" CpGs
            if reInterior + mcrbcInterior == 0:
                continue
            # Rebuild the CpG object
            cpg = CpGFull(coord)
            cpg.re_idx = re_idx
            cpg.mcrbc_idx = mcrbc_idx
            cpg.re_count = re_count
            cpg.mcrbc_count = mcrbc_count
            cpg.re_clean_count = re_clean_count
            cpg.mcrbc_clean_count = mcrbc_clean_count
            cpg.reInterior = reInterior
            cpg.reEnd = reEnd
            cpg.mcrbcInterior = mcrbcInterior
            try:
                self.cpg_sites.add(cpg)
            except MethError, e:
                print e.value
                sys.exit(2)
        cpgfh.close()
                    
    def estimate_methprob(self, p_bar):
        "Estimate methylation probability score (p) and coverage (n) for each CpG."
        if p_bar < EPSILON or abs(p_bar - 1.0) < EPSILON:
            raise MethError('Global methylation level should be between 0 and 1 (0 < meth < 1)')
        sum_n1, sum_n2 = self._sumn_in_segments()
        try:
            p1_p2 = self._get_sampling_paras(p_bar, sum_n1, sum_n2)
        except MethError, e:
            raise MethError(e.value)
        for cpg in self.cpg_sites.sorted_iter():
            try:
                cpg.cal_prob(p1_p2)
            except MethError, e:
                raise MethError(e.value)
        p_bar_new = self._cal_p_bar()
        return p1_p2, p_bar_new
