#!/usr/bin/python
import os
import sys
from optparse import OptionParser
######################################################
## Options and arguments
######################################################
## Usage
usage = "usage: %prog [options] FASTA"
## Version
parser = OptionParser(usage=usage, version="BarcodeCluster 1.0.0")
## Personalized parameters
'''barcode_fasta, cdhit, RAM, threads, outdir, name, threshold'''
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
(options, args) = parser.parse_args()
cdhit = options.CDHIT
RAM = options.RAM
threads = options.Threads
threshold = options.Threshold
outdir = options.Outputdir
name = options.Name
#########################################################
## Checking input fastq file
#########################################################
if not len(args) == 1:
    print('Wrong number of barcode fasta file, only one is allowed')
    sys.exit()
barcode_fasta = args[0]
############################################################
## Function definition
############################################################
## Cluster by CDHIT
def CD_HIT(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold):
    return_number = os.system('%s 1>/dev/null 2>/dev/null' % (cdhit))
    print(return_number)
    if not return_number == 256:
        raise SystemError('CD-HIT cannot be called: %s' % (cdhit))
    else:
        out_fasta = os.path.join(outdir, name+'group.fasta')
        os.system('%s -d 150 -i %s -o %s -c %s -A 1.0 -g 1 -M %s -T %s' % (cdhit, barcode_fasta, out_fasta, threshold, RAM, threads))
        return out_fasta+'.clstr'
## Read the cluster group
def ClusterRead(group_txt):
    clusters = {}
    group_txt_in = list(open(group_txt, 'r'))
    n = 0
    for line in group_txt_in:
        if line.strip().startswith('>'):
            cluster_order = line.strip().strip('>').replace(' ', '_').replace('Cluster', 'Cell')
            clusters[cluster_order] = []
            break_boolean, z, af_lines =  False, 0, group_txt_in[n+1:]
            while break_boolean == False:
                if z < len(af_lines):
                    linex = af_lines[z]
                    if not linex.strip().startswith('>'):
                        tag = linex.strip().split(',')[1].strip().strip('>').split('...')[0]
                        clusters[cluster_order].append(tag)
                        z += 1
                    else:
                        break_boolean = True
                else:
                    break_boolean = True
        n += 1
    return clusters
## Main function
def main(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold):
    group_txt = CD_HIT(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold)
    clusters = ClusterRead(group_txt)
    with open(os.path.join(outdir, name+'.cells_reads.txt'), 'w') as cluster_txt:
        for k,v in clusters.items():
            outline = k+'\t'
            for read in v:
                outline += read + '\t'
            cluster_txt.write(outline.strip() + '\n')
############################################################
## Work flow
############################################################
if (not cdhit == None) and (not outdir == None) and (not name == None):
    main(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold)
else:
    print('Reiquired parameters missing')
    sys.exit()