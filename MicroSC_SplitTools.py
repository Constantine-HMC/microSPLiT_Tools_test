#!/usr/bin/python
import os
import sys
import copy
import linecache
from optparse import OptionParser
######################################################
## Options and arguments
######################################################
## Global variables
src_dir = os.path.split(os.path.realpath(__file__))[0]
## Usage
usage = "usage: %prog [options] FASTQ1 FASTQ2"
## Version
parser = OptionParser(usage=usage, version="MicroSC_SplitTools 1.0.0")
## Personalized parameters
parser.add_option('-s', '--BarcodeModuleSequence', default=None, type="string", action="store", dest='BarcodeModuleSequence',
                    help='Sequence of "Barcode Module", use "N" for region of variable region,'
                         ' for example, NNNNNNNNNNNNGTGACTGATGAANNNNNNCGTAGTTGCAAANNNNNN,'
                         ' deafult = None')
parser.add_option('-b', '--BarcodeModuleCigar', default=None, type="string", action="store", dest='BarcodeModuleCigar',
                    help='Cigar string of Barcode Module, '
                         'use "N" for region of variable region and "T" for fixed technical sequence, '
                         'for example "6B6B12T6B12T6B", deafult = None')
parser.add_option('-g', '--SearchRange', default=0.95, type="float", action="store", dest='SearchRange',
                    help='Fold of presupposed search range for fixed technical sequence'
                         '(based on the length of ecah fixed technical sequence, '
                         'searching will be conducted upstream and downstream "Fold"*"length of "fixed technical '
                         'sequence" bp of the preset position of the fixed technical sequence),  deafult = 1.25')
parser.add_option('-l', '--BarcodePosition', default=None, type="string", action="store", dest='BarcodePosition',
                    help='The position orders of barcode sequence in "Barcode Module", '
                         'for example "2,4,6", default=None')
parser.add_option('-f', '--ReamoveOrder', default=None, type="string", action="store", dest='ReamoveOrder',
                    help='Remove specified barcode region from sequencing reads, use the position order '
                         'to specified barcode region, froma "start" to "end",'
                         'for example "2,6", deafult = None')
parser.add_option('-c', '--CDHIT', default=None, type="string", action="store", dest='CDHIT',
                    help='The path of CD-HIT, deafult = None')
parser.add_option('-r', '--RAM', default=800, type="int", action="store", dest='RAM',
                    help='random access memory, deafult = 800(M)')
parser.add_option('-p', '--Threads', default=1, type="int", action="store", dest='Threads',
                    help='Number of cores to use, deafult = 1')
parser.add_option('-t', '--Threshold', default=0.95, type="float", action="store", dest='Threshold',
                    help='Alignment threshold of clustering barcode sequence, '
                         'Please calculate according to the allowable mismatch number, deafult = 0.95')
parser.add_option('-o', '--Outputdir', default=None, type="string", action="store", dest='Outputdir',
                    help='Absolute path of outputdir, only reads name and barcode, deafult = None')
parser.add_option('-n', '--Name', default=None, type="string", action="store", dest='Name',
                    help='Prefix name of output files, deafult = None')
## Parse options
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
cdhit = options.CDHIT
RAM = options.RAM
threads = options.Threads
#########################################################
## Checking input fastq file
#########################################################
if not len(args) == 2:
    print('Wrong number of barcode fasta file, only one is allowed')
    sys.exit()
fq1, fq2 = args[0], args[1]
#########################################################
## Other scripts
#########################################################
BarcodeCluster = os.path.join(src_dir, 'BarcodeCluster.py')
ExtractBarcode = os.path.join(src_dir, 'ExtractBarcode.py')
############################################################
## Function definition
############################################################
def BarcoeFastaReadIn(barcode_faa, direction):
    tag_barcode = {}
    with open(barcode_faa, 'r') as in_faa:
        for line in in_faa:
            if line.strip().startswith('>'):
                tag = line.strip().strip('>')
                barcode = in_faa.readline().strip()
                tag_barcode[tag+'_lzc_'+direction] = barcode
    return tag_barcode
def BarcoeFastaRedup(barcode_faa1, barcode_faa2):
    tag_barcode1 = BarcoeFastaReadIn(barcode_faa1, '1')
    tag_barcode2 = BarcoeFastaReadIn(barcode_faa2, '2')
    tag_barcode = copy.deepcopy(tag_barcode1)
    tag_barcode.update(tag_barcode2)
    return tag_barcode
def FastqReadIn(fastq):
    fastq_read = {}
    with open(fastq, 'r') as in_fastq:
        n = 0
        for line in in_fastq:
            n += 1
            if line.strip().startswith('@'):
                if not linecache.getline(fastq, n-1).strip().startswith('+'):
                    seq = linecache.getline(fastq, n + 1).strip()
                    qual = linecache.getline(fastq, n + 3).strip()
                    seq_tag = line.strip().strip('@').split(' ')[0]
                    fastq_read[seq_tag] = [seq, qual]
    return fastq_read
############################################################
## Work flow
############################################################
# print('python3 %s -s %s -c %s -t %s -r %s -p %s -f %s -o %s -n %s %s' % (ExtractBarcode, module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, outdir, name+'_1', fq1))
# print('python3 %s -s %s -c %s -t %s -r %s -p %s -f %s -o %s -n %s %s' % (ExtractBarcode, module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, outdir, name+'_2', fq2))
# os.system('python3 %s -s %s -c %s -t %s -r %s -p %s -f %s -o %s -n %s %s' % (ExtractBarcode, module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, outdir, name+'_1', fq1))
barcode_fq1 = os.path.join(outdir, name+'_1' + '.RemoveBarcode.fastq')
barcode_faa1 = os.path.join(outdir, name+'_1' + '.Barcode.fasta')
# os.system('python3 %s -s %s -c %s -t %s -r %s -p %s -f %s -o %s -n %s %s' % (ExtractBarcode, module_seq, module_cigar, threshold, range_fold, barcode_position, remove_order, outdir, name+'_2', fq2))
barcode_fq2 = os.path.join(outdir, name+'_2' + '.RemoveBarcode.fastq')
barcode_faa2 = os.path.join(outdir, name+'_2' + '.Barcode.fasta')
tag_barcode = BarcoeFastaRedup(barcode_faa1, barcode_faa2)
PE_barcode_faa = os.path.join(outdir, name+'.PE.Barcode.fasta')
with open(PE_barcode_faa, 'w') as out1:
    for tag,barcode in tag_barcode.items():
        out1.write('>%s\n%s\n' % (tag, barcode))
os.system('python3 %s -c %s -r %s -p %s -t %s -o %s -n %s %s' % (BarcodeCluster, cdhit, RAM, threads, threshold, outdir, name, PE_barcode_faa))
cell_reads_txt = os.path.join(outdir, name+'.cells_reads.txt')
cells_reads = {}
with open(cell_reads_txt, 'r') as in1:
    for line in in1:
        ll = line.strip().split('\t')
        cell, reads = ll[0], ll[1:]
        cells_reads[cell] = reads
row_fastq1_read = FastqReadIn(fq1)
print(len(row_fastq1_read))
row_fastq2_read = FastqReadIn(fq2)
print(len(row_fastq2_read))
barcode_fastq1_read = FastqReadIn(barcode_fq1)
print(len(barcode_fastq1_read))
barcode_fastq2_read = FastqReadIn(barcode_fq2)
print(len(barcode_fastq2_read))
fq_dir = os.path.join(outdir, 'cell_fastq')
os.mkdir(fq_dir)
for cell,reads in cells_reads.items():
    cell_dir = os.path.join(fq_dir, cell)
    os.mkdir(cell_dir)
    cell_fq1 = open(os.path.join(cell_dir, '1.fastq'), 'w')
    cell_fq2 = open(os.path.join(cell_dir, '2.fastq'), 'w')
    read_pair = {}
    for tag_source in reads:
        tag = tag_source.strip().split('_lzc_')[0]
        source = tag_source.strip().split('_lzc_')[1]
        if not tag in read_pair.keys():
            read_pair[tag] = 'IN'
            if source == '1':
                cell_fq1.write('@%s\n%s\n+\n%s\n' % (tag, barcode_fastq1_read[tag][0], barcode_fastq1_read[tag][1]))
                if tag in barcode_fastq2_read.keys():
                    cell_fq2.write('@%s\n%s\n+\n%s\n' % (tag, barcode_fastq2_read[tag][0], barcode_fastq2_read[tag][1]))
                elif not tag in barcode_fastq2_read.keys() and tag in row_fastq2_read.keys():
                    cell_fq2.write('@%s\n%s\n+\n%s\n' % (tag, row_fastq2_read[tag][0], row_fastq2_read[tag][1]))
            elif source == '2':
                cell_fq2.write('@%s\n%s\n+\n%s\n' % (tag, barcode_fastq2_read[tag][0], barcode_fastq2_read[tag][1]))
                if tag in barcode_fastq1_read.keys():
                    cell_fq1.write('@%s\n%s\n+\n%s\n' % (tag, barcode_fastq1_read[tag][0], barcode_fastq1_read[tag][1]))
                elif not tag in barcode_fastq1_read.keys() and tag in row_fastq1_read.keys():
                    cell_fq1.write('@%s\n%s\n+\n%s\n' % (tag, row_fastq1_read[tag][0], row_fastq1_read[tag][1]))
