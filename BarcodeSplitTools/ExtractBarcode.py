#!/usr/bin/python
import re
import os
import sys
from fuzzywuzzy import fuzz
import linecache
from optparse import OptionParser

######################################################
## 选项和参数
######################################################
## 用法
usage = "Usage: %prog [options] FASTQ\n" \
        "The extracted barcode sequence will be filled by A if barcode length less than 18" \
        "\nBarcode sequence will be removed from sequencing reads (the maximum coverage of the input position orders of barcode and anchor)"
## 版本
parser = OptionParser(usage=usage, version='ExtractBarcode.py test_version:0.0.1')
## 个性化选项
parser.add_option('-m', '--module', default=None, type="string", action="store", dest="module", metavar='string',
                    help='Required, sequence of "Barcode Module", use "N" for variable region, for example, "NNNNNNNNNNNNGTGACTGATGAANNNNNNCGTAGTTGCAAANNNNNN", deafult = None')
parser.add_option('-c', '--cigar', default=None, type="string", action="store", dest='cigar', metavar='string',
                    help='Required, cigar string of "Barcode Module", use "N" variable region and "T" for fixed anchoring technical sequence, for example "6B6B12T6B12T6B", deafult = None')
parser.add_option('-b', '--barcode', default=None, type="string", action="store", dest='barcode', metavar='string',
                    help='Required, the position orders of barcode sequence in "ModuleCigar", for example "2,4,6", default = None')
parser.add_option('-s', '--search', default='False', type="string", action="store", dest='search', metavar='boolean',
                    help='Optional, to determine the position of barcode by fixed anchoring technical sequence, enter "True" or "False", default = False')
parser.add_option('-a', '--anchor', default=None, type="string", action="store", dest='anchor', metavar='string',
                    help='Optional but requierd when "search" if used, should be used with "-s", the position orders of anchoring technical sequence in "Barcode Module", for example "3,5", default = None')
parser.add_option('-f', '--fold', default=0.95, type="float", action="store", dest='fold', metavar='float',
                    help='Optional, should be used with "-s", fold of presupposed search range for fixed anchoring technical sequence, deafult = 1.25')
parser.add_option('-t', '--threshold', default=0.95, type="float", action="store", dest='threshold', metavar='float',
                    help='Optional, should be used with "-s", alignment threshold of fixed anchoring technical sequence, deafult = 0.95')
parser.add_option('-o', '--outputdir', default=None, type="string", action="store", dest='outputdir', metavar='path',
                    help='Required, absolute path of outputdir, only reads name and barcode, deafult = None')
parser.add_option('-n', '--name', default=None, type="string", action="store", dest='name', metavar='string',
                    help='Required, prefix name of output files, deafult = None')
## 解析输入的选项
(options, args) = parser.parse_args()
module_seq = options.module
module_cigar = options.cigar
search_tf = options.search
threshold = options.threshold
range_fold = options.fold
barcode_position = options.barcode
anchor_position = options.anchor
outdir = options.outputdir
name = options.name
#########################################################
## 检查部分选项的输入
#########################################################
## 检查是否使用锚定序列
if not search_tf in ['True', 'False']:
    print('Wrong boolean of "search", please specify "True" of "False"!')
    sys.exit()
else:
    search_tf = True if search_tf == "True" else False
## 检查使用锚定序列时所需要的其他输入
if search_tf == True:
    if anchor_position == None:
        print('Please specify the position orders of anchoring technical sequence in "Barcode Module"!')
        sys.exit()
#########################################################
## 检查输入的FASTQ文件
#########################################################
if not len(args) == 1:
    print('Wrong number of fastq file, only one fastq file is allowed')
    sys.exit()
fastq = args[0]
############################################################
## 函数定义
############################################################
## 读取模板序列的cigar
def ModuleSeqSplit(module_cigar, module_seq):
    ## 读取cigar中的序列类型
    subseq_type = [x for x in re.split(r'\d+', module_cigar)[1:]]
    ## 读取cigar中的序列长度
    subseq_len = [int(x) for x in re.split(r'[a-zA-Z]', module_cigar)[0:-1]]
    ## 计算每一段序列的起点位置
    accumulated_len = [sum(subseq_len[0:x]) for x in list(range(0,len(subseq_len)))]
    ## 读取每一段序列
    subseq = [module_seq[accumulated_len[x]:accumulated_len[x]+subseq_len[x]] for x in list(range(0, len(subseq_type)))]
    ## 合并每段序列的多个信息
    module_seq_split = [list(x) for x in zip(subseq_type, accumulated_len, subseq_len, subseq)]
    return module_seq_split
## 读取指定的序列位置
def PositionSplit(position):
    ## 切片字符串并转化为列表
    order_list = position.strip().split(',')
    try:
        order_list = [int(x)-1 for x in order_list]
    except:
        print('Wrong input of position order!')
        sys.exit()
    return order_list
## 读取FASTQ文件
def FastqReadIn(fastq):
    txt_len = len(open(fastq, 'r').readlines())
    fastq_read, read_count, n = {}, int(txt_len/4), 1
    while n < txt_len:
        if linecache.getline(fastq, n).strip().startswith('@') and linecache.getline(fastq, n+2).strip().startswith('+'):
            seq_tag = linecache.getline(fastq, n).strip().strip('@').split(' ')[0]
            seq = linecache.getline(fastq, n + 1).strip()
            qual = linecache.getline(fastq, n + 3).strip()
            fastq_read[seq_tag] = [seq, qual]
        n += 1
    return fastq_read, read_count
## 序列转化kmer
def SeqToKmers(seq, kmer_len):
    kmers = []
    for x in list(range(0,len(seq)-kmer_len+1)):
        kmers.append(seq[x:x+kmer_len])
    return kmers
## 需要从测序reads种去除的范围
def RemoveRange(module_seq_split, barcode_order_list, anchor_order_list=None):
    ## 基本的去除范围：barcode序列的所覆盖的范围
    all_remove_region = barcode_order_list
    ## 检查是否输入锚定技术序列，将锚定技术序列加入去除范围
    if not anchor_order_list == None:
        all_remove_region = list(set(barcode_order_list + anchor_order_list))
    remove_start = module_seq_split[min(all_remove_region)][1]
    remove_end = module_seq_split[max(all_remove_region)][1] + module_seq_split[max(all_remove_region)][2]
    return [remove_start, remove_end]
## 使用锚定技术序列锚定的情况下提取Barcode
def BarcodeExtractAnchor(seq, qual, module_seq_split, barcode_order_list, anchor_order_list, threshold, range_fold):
    ## 字典，录入所有锚定技术序列的锚定位置
    tech_anchor_position = {}
    ## 循环，寻找所有锚定序列的锚定位置
    for anchor_tech in  anchor_order_list:
        ## 读取锚定序列的相关信息
        tech_seq_inf = module_seq_split[anchor_tech]
        ## 字典，录入在搜索锚定序列的范围内，锚定序列的可能锚定位置
        align_positions = {}
        ## 根据预设的range_fold确定搜索锚定技术序列的的范围
        search_range = [max(0, int(tech_seq_inf[1] - tech_seq_inf[2] * range_fold)),
                        min(len(seq), int(tech_seq_inf[1] + tech_seq_inf[2] + tech_seq_inf[2] * range_fold))]
        ## 搜索范围的序列打断为kmer
        kmers = SeqToKmers(seq[search_range[0]:search_range[1] + 1], tech_seq_inf[2])
        ## 搜索锚定序列的锚定位置
        for pos_index,kmer in enumerate(kmers):
            align_ratio = fuzz.ratio(kmer, tech_seq_inf[3])
            if align_ratio >= threshold * 100:
                align_position = pos_index + search_range[0] - tech_seq_inf[1]
                align_positions[align_position] = align_positions[align_position] + align_ratio if align_position in align_positions.keys() else align_ratio
        ## 确认锚定序列的锚定位置
        exact_align_position = sorted(align_positions.items(), key=lambda x: x[1], reverse=True)[0][0] if not align_positions == {} else None
        tech_anchor_position[anchor_tech] = exact_align_position
    ## 确认多个锚定技术序列是否使read整体具有相同的比对位置
    if len(list(set(list(tech_anchor_position.values())))) == 1 and not None in tech_anchor_position.values():
        exact_anchor_position = list(tech_anchor_position.values())[0]
        aligned_seq = seq[exact_anchor_position:] if exact_anchor_position >=0 else 'N'*abs(exact_anchor_position)+seq
        ## 获取多个barcode区域序列并组合在一起
        barcode = ''
        for barcode_region in barcode_order_list:
            barcode += aligned_seq[module_seq_split[barcode_region][1]: module_seq_split[barcode_region][1] +module_seq_split[barcode_region][2]]
        if len(barcode) < 18:
            barcode += (18-len(barcode))*'A'
        ## 要去除的范围
        remove_range = RemoveRange(module_seq_split, barcode_order_list, anchor_order_list)
        ## 去除barcode区域，和锚定技术序列
        new_seq, new_qual = None, None
        if exact_anchor_position >= 0 or (exact_anchor_position < 0 and remove_range[0] + exact_anchor_position >= 0):
            new_seq = seq[0:remove_range[0] + exact_anchor_position] + seq[remove_range[1] + exact_anchor_position:]
            new_qual = qual[0:remove_range[0] + exact_anchor_position] + qual[remove_range[1] + exact_anchor_position:]
        elif exact_anchor_position < 0 and remove_range[0] + exact_anchor_position < 0:
            new_seq = seq[remove_range[1] + exact_anchor_position:]
            new_qual = qual[remove_range[1] + exact_anchor_position:]
        if (not new_seq == None) and (not new_qual == None) and (not barcode == ''):
            return barcode, new_seq, new_qual
        else:
            return None
    else:
        return None
## 不使用锚定技术序列锚定的情况下提取Barcode
def BarcodeExtractNoneAnchor(seq, qual, module_seq_split, barcode_order_list, anchor_order_list):
    ## 获取多个barcode区域序列并组合在一起
    barcode = ''
    for barcode_region in barcode_order_list:
        barcode += seq[module_seq_split[barcode_region][1]: module_seq_split[barcode_region][1] + module_seq_split[barcode_region][2]]
    if len(barcode) < 18:
        barcode += (18 - len(barcode)) * 'A'
    ## 去除范围
    remove_range = RemoveRange(module_seq_split, barcode_order_list, anchor_order_list)
    ## 去除barcode区域，和锚定技术序列
    exact_anchor_position = 0
    new_seq, new_qual = None, None
    if exact_anchor_position >= 0 or (exact_anchor_position < 0 and remove_range[0] + exact_anchor_position >= 0):
        new_seq = seq[0:remove_range[0] + exact_anchor_position] + seq[remove_range[1] + exact_anchor_position:]
        new_qual = qual[0:remove_range[0] + exact_anchor_position] + qual[remove_range[1] + exact_anchor_position:]
    elif exact_anchor_position < 0 and remove_range[0] + exact_anchor_position < 0:
        new_seq = seq[remove_range[1] + exact_anchor_position:]
        new_qual = qual[remove_range[1] + exact_anchor_position:]
    if (not new_seq == None) and (not new_qual == None) and (not barcode == ''):
        return barcode, new_seq, new_qual
    else:
        return None
## 主函数
def main(module_seq, module_cigar, threshold, range_fold, barcode_position, anchor_position, fastq, outdir, name, search=True):
    ## 读入fastq文件
    fastq_read, read_count = FastqReadIn(fastq)
    ## 读取barcode序列模式
    module_seq_split = ModuleSeqSplit(module_cigar, module_seq)
    ## 读取barcode区域
    barcode_order_list = PositionSplit(barcode_position)
    ## 读取技术序列区域
    anchor_order_list = None
    if not anchor_position == None:
        anchor_order_list = PositionSplit(anchor_position)
    ## 输出文件
    barcode_fasta = open(os.path.join(outdir, name + '.Barcode.fasta'), 'w')
    remove_barcode_fastq = open(os.path.join(outdir, name + '.RemoveBarcode.fastq'), 'w')
    no_barcode = open(os.path.join(outdir, name + '.Reads.NonBarcoded.txt'), 'w')
    summary = open(os.path.join(outdir, name + '.ExtractBarcode.summary.txt'), 'w')
    ## 无标签的序列数量
    nonbar_count = 0
    ## 循环， 处理每一个read
    for seq_tag, inf in fastq_read.items():
        ## 具体的序列，碱基质量编码
        seq, qual = inf[0], inf[1]
        ## 进行锚定序列锚定，提取barcode并去除barcode区域，和锚定技术序列
        if search == True:
            extract_results = BarcodeExtractAnchor(seq, qual, module_seq_split, barcode_order_list, anchor_order_list, threshold, range_fold)
            if extract_results != None:
                barcode, new_seq, new_qual = extract_results[0], extract_results[1], extract_results[2]
                barcode_fasta.write('>%s\n%s\n' % (seq_tag, barcode))
                remove_barcode_fastq.write('@%s\n%s\n+\n%s\n' % (seq_tag, new_seq, new_qual))
            elif extract_results == None:
                nonbar_count += 1
                no_barcode.write(seq_tag + '\n')
        ## 不进行锚定序列锚定，提取barcode并去除barcode区域，和锚定技术序列
        elif search == False:
            extract_results = BarcodeExtractNoneAnchor(seq, qual, module_seq_split, barcode_order_list, anchor_order_list)
            if extract_results != None:
                barcode, new_seq, new_qual = extract_results[0], extract_results[1], extract_results[2]
                barcode_fasta.write('>%s\n%s\n' % (seq_tag, barcode))
                remove_barcode_fastq.write('@%s\n%s\n+\n%s\n' % (seq_tag, new_seq, new_qual))
            elif extract_results == None:
                nonbar_count += 1
                no_barcode.write(seq_tag + '\n')
    ## 输出summary.txt
    summary.write('%s\nTotal Reads: %s\nNonBarcoded Reads Count: %s\nNonBarcoded Reads Ratio: %s' % (fastq, read_count, nonbar_count, str(round(nonbar_count / read_count * 100, 4)) + '%'))
############################################################
## 工作流程
############################################################
if (not module_seq == None) and (not module_cigar == None) and (not barcode_position == None) and (not fastq == None) and (not outdir == None) and (not name == None):
    main(module_seq, module_cigar, threshold, range_fold, barcode_position, anchor_position, fastq, outdir, name, search_tf)
else:
    print('Reiquired parameters missing')
    sys.exit()