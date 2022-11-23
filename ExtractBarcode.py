#!/usr/bin/python
import re
import os
import sys
from fuzzywuzzy import fuzz
import linecache
from optparse import OptionParser

######################################################
## Options and arguments
######################################################
## Usage
usage = "usage: %prog [options] FASTQ"
## Version
parser = OptionParser(usage=usage, version="ExtractBarcode 1.0.0")
## Personalized parameters
parser.add_option('-s', '--BarcodeModuleSequence', default=None, type="string", action="store", dest='BarcodeModuleSequence',
                    help='Sequence of "Barcode Module", use "N" for region of variable region,'
                         ' for example, NNNNNNNNNNNNGTGACTGATGAANNNNNNCGTAGTTGCAAANNNNNN,'
                         ' deafult = None')
parser.add_option('-c', '--BarcodeModuleCigar', default=None, type="string", action="store", dest='BarcodeModuleCigar',
                    help='Cigar string of Barcode Module, '
                         'use "N" for region of variable region and "T" for fixed technical sequence, '
                         'for example "6B6B12T6B12T6B", deafult = None')
parser.add_option('-t', '--Threshold', default=0.95, type="float", action="store", dest='Threshold',
                    help='Alignment threshold of fixed technical sequence, deafult = 0.95')
parser.add_option('-r', '--SearchRange', default=0.95, type="float", action="store", dest='SearchRange',
                    help='Fold of presupposed search range for fixed technical sequence'
                         '(based on the length of ecah fixed technical sequence, '
                         'searching will be conducted upstream and downstream "Fold"*"length of "fixed technical '
                         'sequence" bp of the preset position of the fixed technical sequence),  deafult = 1.25')
parser.add_option('-p', '--BarcodePosition', default=None, type="string", action="store", dest='BarcodePosition',
                    help='The position orders of barcode sequence in "Barcode Module", '
                         'for example "2,4,6", default=None')
parser.add_option('-f', '--ReamoveOrder', default=None, type="string", action="store", dest='ReamoveOrder',
                    help='Remove specified barcode region from sequencing reads, use the position order '
                         'to specified barcode region, froma "start" to "end",'
                         'for example "2,6", deafult = None')
parser.add_option('-o', '--Outputdir', default=None, type="string", action="store", dest='Outputdir',
                    help='Absolute path of outputdir, only reads name and barcode, deafult = None')
parser.add_option('-n', '--Name', default=None, type="string", action="store", dest='Name',
                    help='Prefix name of output files, deafult = None')
## Parse options
(options, args) = parser.parse_args()
module_seq = options.BarcodeModuleSequence
module_cigar = options.BarcodeModuleCigar
threshold = options.Threshold
range_fold = options.SearchRange
barcode_position = options.BarcodePosition
remove_order = options.ReamoveOrder
outdir = options.Outputdir
name = options.Name
#########################################################
## Checking input fastq file
#########################################################
if not len(args) == 1:
    print('Wrong number of fastq file, only one fastq file is allowed')
    sys.exit()
fastq = args[0]
############################################################
## Function definition
############################################################
## Split barcode sequence module
def ModuleSeqSplit(module_cigar, module_seq):
    subseq_type = [x for x in re.split(r'\d+', module_cigar)[1:]]
    subseq_len = [int(x) for x in re.split(r'[a-zA-Z]', module_cigar)[0:-1]]
    accumulated_len = [sum(subseq_len[0:x]) for x in list(range(0,len(subseq_len)))]
    subseq = [module_seq[accumulated_len[x]:accumulated_len[x]+subseq_len[x]] for x in list(range(0, len(subseq_type)))]
    module_seq_split = [list(x) for x in zip(subseq_type, accumulated_len, subseq_len, subseq)]
    return module_seq_split
## Split specified barcode position
def BarcodePositionSplit(barcode_position):
    barcode_order_list = barcode_position.strip().split(',')
    barcode_order_list = [int(x)-1 for x in barcode_order_list]
    return barcode_order_list
## Read fastq file
def FastqReadIn(fastq):
    fastq_read = {}
    read_count = 0
    with open(fastq, 'r') as in_fastq:
        n = 0
        for line in in_fastq:
            n += 1
            if line.strip().startswith('@'):
                if not linecache.getline(fastq, n-1).strip().startswith('+'):
                    read_count += 1
                    seq = linecache.getline(fastq, n + 1).strip()
                    qual = linecache.getline(fastq, n + 3).strip()
                    seq_tag = line.strip().strip('@').split(' ')[0]
                    fastq_read[seq_tag] = [seq, qual]
    return fastq_read, read_count
## Split sequence to kmers
def SeqToKmers(seq, kmer_len):
    kmers = []
    for x in list(range(0,len(seq)-kmer_len+1)):
        kmers.append(seq[x:x+kmer_len])
    return kmers
## Extract barcode sequence from sequencing read
def BarcodeExtract(test_seq, module_seq_split, barcode_order_list, threshold, range_fold):
    tech_seq_pos = []
    tech_seq_check = {}
    tech_count = 0
    for seq_region in module_seq_split:
        if seq_region[0] == 'T':
            tech_count += 1
            tech_seq_check[tech_count] = False
            align_positions = {}
            search_range = [max(0, int(seq_region[1]-seq_region[2]*range_fold)), min(len(test_seq), int(seq_region[1]+seq_region[2]+seq_region[2]*range_fold))]
            kmers = SeqToKmers(test_seq[search_range[0]:search_range[1]+1], seq_region[2])
            for index,kmer in enumerate(kmers):
                align_ratio = fuzz.ratio(kmer, seq_region[3])
                if align_ratio >= threshold*100:
                    tech_seq_check[tech_count] = True
                    align_position = index + search_range[0] - seq_region[1]
                    if align_position in align_positions:
                        align_positions[align_position] += align_ratio
                    elif not align_position in align_positions:
                        align_positions[align_position] = align_ratio
            if not align_positions == {}:
                exact_align_position = sorted(align_positions.items(), key=lambda x: x[1], reverse=True)[0][0]
                if not exact_align_position in tech_seq_pos:
                    tech_seq_pos.append(exact_align_position)
    if len(tech_seq_pos) == 1:
        if not False in tech_seq_check.values():
            if tech_seq_pos[0] >= 0:
                align_seq = test_seq[tech_seq_pos[0]:]
            else:
                align_seq = 'N'*abs(tech_seq_pos[0])+test_seq
            barcode = ''
            for barcode_oder in barcode_order_list:
                barcode += align_seq[module_seq_split[barcode_oder][1]: module_seq_split[barcode_oder][1]+module_seq_split[barcode_oder][2]]
            return barcode, tech_seq_pos
        else:
            return None
    else:
        return None
## Get the barcode and technical sequence region that should be removed
def RemoveRange(remove_order, module_seq_split):
    remove_order_list = remove_order.strip().split(',')
    remove_order_list = [int(x) - 1 for x in remove_order_list]
    remove_start = module_seq_split[remove_order_list[0]][1]
    remove_end = module_seq_split[remove_order_list[1]][1] + module_seq_split[remove_order_list[1]][2]
    return [remove_start, remove_end]
## main funtion
def main(module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, fastq, outdir, name):
    fastq_read, read_count = FastqReadIn(fastq)
    module_seq_split = ModuleSeqSplit(module_cigar, module_seq)
    remove_range = RemoveRange(remove_order, module_seq_split)
    barcode_order_list = BarcodePositionSplit(barcode_position)
    barcode_fasta = open(os.path.join(outdir, name + '.Barcode.fasta'), 'w')
    no_barcode = open(os.path.join(outdir, name + '.Reads.NonBarcoded.txt'), 'w')
    remove_barcode_fastq = open(os.path.join(outdir, name + '.RemoveBarcode.fastq'), 'w')
    nonbar_count = 0
    for seq_tag, inf in fastq_read.items():
        seq, qual = inf[0], inf[1]
        barcode_align_pos = BarcodeExtract(seq, module_seq_split, barcode_order_list, threshold, range_fold)
        if not barcode_align_pos == None:
            barcode = barcode_align_pos[0]
            tech_seq_pos = barcode_align_pos[1][0]
            if tech_seq_pos >= 0:
                new_seq = seq[0:remove_range[0] + tech_seq_pos] + seq[remove_range[1] + tech_seq_pos:]
                new_qual = qual[0:remove_range[0] + tech_seq_pos] + qual[remove_range[1] + tech_seq_pos:]
                remove_barcode_fastq.write('@%s\n%s\n+\n%s\n' % (seq_tag, new_seq, new_qual))
            elif tech_seq_pos < 0:
                if remove_range[0] + tech_seq_pos < 0:
                    new_seq = seq[remove_range[1] + tech_seq_pos:]
                    new_qual = qual[remove_range[1] + tech_seq_pos:]
                    remove_barcode_fastq.write('@%s\n%s\n+\n%s\n' % (seq_tag, new_seq, new_qual))
                elif remove_range[0] + tech_seq_pos >= 0:
                    new_seq = seq[0:remove_range[0] + tech_seq_pos] + seq[remove_range[1] + tech_seq_pos:]
                    new_qual = qual[0:remove_range[0] + tech_seq_pos] + qual[remove_range[1] + tech_seq_pos:]
                    remove_barcode_fastq.write('@%s\n%s\n+\n%s\n' % (seq_tag, new_seq, new_qual))
            barcode_fasta.write('>%s\n%s\n' % (seq_tag, barcode))
        else:
            nonbar_count += 1
            no_barcode.write(seq_tag+'\n')
    print('%s\nTotal Reads: %s\nNonBarcoded Reads Count: %s\nNonBarcoded Reads Ratio: %s' % (fastq, read_count, nonbar_count, str(round(nonbar_count/read_count*100, 4))+'%'))
############################################################
## Work flow
############################################################
if (not module_seq == None) and (not module_cigar == None) and (not barcode_position == None) and (not fastq == None) and (not outdir == None):
    main(module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, fastq, outdir, name)
else:
    print('Reiquired parameters missing')
    sys.exit()