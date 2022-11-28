#!/usr/bin/python
import os
import sys
from optparse import OptionParser
######################################################
## 选项和参数
######################################################
## 用法
usage = "usage: %prog [options] FASTA"
## 版本
parser = OptionParser(usage=usage, version="BarcodeCluster test_version:0.0.1")
## 个性化选项
parser.add_option('-c', '--CD-HIT', default=None, type="string", action="store", dest='cdhit', metavar='path',
                    help='Required, The path of cd-hit-est, deafult = None')
parser.add_option('-r', '--RAM', default=800, type="int", action="store", dest='RAM', metavar='integer',
                    help='Optional, random access memory, deafult = 800(M)')
parser.add_option('-p', '--threads', default=1, type="int", action="store", dest='threads', metavar='integer',
                    help='Optional, Number of threads to use, deafult = 1')
parser.add_option('-t', '--threshold', default=0.95, type="float", action="store", dest='threshold', metavar='float',
                    help='Optional, Alignment threshold of clustering barcode sequence, Please calculate according to the allowable mismatch number, deafult = 0.95')
parser.add_option('-o', '--outputdir', default=None, type="string", action="store", dest='outputdir', metavar='path',
                    help='Required, Absolute path of outputdir, only reads name and barcode, deafult = None')
parser.add_option('-n', '--name', default=None, type="string", action="store", dest='name', metavar='string',
                    help='Required, Prefix name of output files, deafult = None')
## 解析输入的选项
(options, args) = parser.parse_args()
cdhit = options.cdhit
RAM = options.RAM
threads = options.threads
threshold = options.threshold
outdir = options.outputdir
name = options.name
#########################################################
## 检查输入的FASTA文件
#########################################################
if not len(args) == 1:
    print('Wrong number of barcode fasta file, only one is allowed')
    sys.exit()
barcode_fasta = args[0]
############################################################
## 函数定义
############################################################
## 用CD-HIT对barcode进行分类，确定细胞数目
def CD_HIT(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold):
    return_number = os.system('%s 1>/dev/null 2>/dev/null' % (cdhit))
    if not return_number == 256:
        raise SystemError('CD-HIT cannot be called: %s' % (cdhit))
    else:
        out_fasta = os.path.join(outdir, name+'.group.fasta')
        os.system('%s -d 150 -i %s -o %s -c %s -A 1.0 -g 1 -M %s -T %s 1>/dev/null 2>/dev/null' % (cdhit, barcode_fasta, out_fasta, threshold, RAM, threads))
        return out_fasta+'.clstr'
## 获取每个分类（每个细胞）下的reads标签
def ClusterRead(group_txt):
    clusters = {}
    group_txt_in, n = list(open(group_txt, 'r')), 0
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
## 主函数
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
## 工作流程
############################################################
if (not cdhit == None) and (not outdir == None) and (not name == None):
    main(barcode_fasta, cdhit, RAM, threads, outdir, name, threshold)
else:
    print('Reiquired parameters missing')
    sys.exit()