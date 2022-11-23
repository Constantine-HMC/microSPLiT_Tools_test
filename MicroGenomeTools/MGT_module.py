"""MGT_module(MicroGenomeTools_module)"""
import shutil

"""
authors:
Zhuochong Liu, Flavius Constantine, Motochika Hata

version:
modudel_test.0.0.1

description:
    This is MGT_module, all function modules of MicroGenomeTools are here.
"""

import os

"""基础函数模块"""
## 检查文件是否存在
def check_file(file):
    abs_path=os.path.abspath(file)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError('No such file: %s' % (abs_path))

## 检查文件夹是否存在，并且创建文件夹
def check_dir(dir, create, clean):
    abs_path=os.path.abspath(dir).rstrip('/')
    if not os.path.isdir(abs_path):
        if create==False:
            raise FileExistsError('Dir not exists: %s' % (abs_path))
        elif create==True:
            if not os.path.isdir(os.path.dirname(abs_path)):
                raise FileExistsError('Dir not exists: %s' % (os.path.dirname(abs_path)))
            else:
                os.mkdir(abs_path)
    elif os.path.isdir(abs_path):
        if clean==True:
            shutil.rmtree(abs_path)
            os.mkdir(abs_path)

## 检查分析流程涉及的软件是否可调用
def check_software(software_list):
    for software in software_list:
        return_number = os.system('%s 1>/dev/null 2>/dev/null' % (software))
        if not (return_number == 0 or return_number == 256 or return_number == 65280):
            print(software, return_number)
            raise SystemError('Software cannot be called: %s' % (software))

## 为软件运行的命令行添加参数
def join_parameter(cmd_head, parameters):
    parameter_str = ''
    for parameter in parameters:
        pl = parameters[parameter].strip().split(':')
        if pl[0].strip() == 's':
            parameter_str += '-' + parameter + ' ' + pl[1].strip().replace('TF','') + ' '
        elif pl[0].strip() == 'a':
            parameter_str += '--' + parameter + ' ' + pl[1].strip().replace('TF','') + ' '
        elif pl[0].strip() == 'e':
            parameter_str += parameter + '=' + pl[1].strip().replace('TF','') + ' '
        else:
            raise TypeError('Wrong parameter type: %s: %s' % (cmd_head, parameter))
    return parameter_str

## 清除中间文件
def clean(file_list, shell_script):
    clean_cmd='rm -fr'
    for file in file_list:
        clean_cmd+=' '+file
    shell_script += clean_cmd+'\n'
    return shell_script

"""流程通用步骤函数模块"""
## fastqc质控
def fq_fastqc(myconf, fastqc, fastq, outdir, threads, shell_scripts):
    fastqc_cmd = '%s ' % (fastqc)
    if 'fastqc' in myconf:
        fastqc_cmd += join_parameter('fastqc', myconf['fastqc'])
    fastqc_cmd += '-t %s -o %s %s %s\n' % (threads, outdir, fastq[0], fastq[1])
    shell_scripts += fastqc_cmd
    return shell_scripts

## fastp质控
def fq_fastp(myconf, fastp, name, fastq, outdir, threads, shell_script):
    abs_prefix1 = os.path.basename(os.path.abspath(fastq[0])).strip('fastq').strip('fq').strip('fastq.gz').strip('fq.gz')
    abs_prefix2 = os.path.basename(os.path.abspath(fastq[1])).strip('fastq').strip('fq').strip('fastq.gz').strip('fq.gz')
    fastp_cmd = '%s -w %s ' % (fastp, threads)
    if 'fastp' in myconf:
        fastp_cmd += join_parameter('fastp', myconf['fastp'])
    qcfq1, qcfq2 = os.path.join(outdir, abs_prefix1 + '.qc.fastq'), os.path.join(outdir, abs_prefix2 + '.qc.fastq')
    html, json = os.path.join(outdir, name + '.fastp.html'), os.path.join(outdir, name + '.fastp.json')
    fastp_cmd += '-j %s -h %s -I %s -O %s -i %s -o %s\n' % (json, html, fastq[0], qcfq1, fastq[1], qcfq2)
    shell_script += fastp_cmd
    return shell_script, [qcfq1, qcfq2]

##  检查参考基因组的系列索引文件是否存在（bwa2）
def ref_index(shell_script, ref, bwa2, samtools):
    bwa2indexfile = [ref + '.0123', ref + '.amb', ref + '.ann', ref + '.bwt.2bit.64', ref + '.pac']
    samtoolsfaidxfile = ref + '.fai'
    for bwa2file in bwa2indexfile:
        if not os.path.isfile(bwa2file):
            shell_script += '%s index %s \n' % (bwa2, ref)
            break
    if not os.path.isfile(samtoolsfaidxfile):
        shell_script += '%s faidx %s \n' % (samtools, ref)
    return shell_script

##  检查参考基因组的系列字典文件是否存在
def ref_dict(shell_script, ref, gatk):
    gatkdict = os.path.join(os.path.dirname(ref), os.path.basename(ref).split('.')[0]) + '.dict'
    if not os.path.isfile(gatkdict):
        shell_script += '%s CreateSequenceDictionary -R %s \n' % (gatk, ref)
    return shell_script

##  检查参考基因组的系列索引文件是否存在（STAR）
def ref_index_STAR(shell_script, myconf, star_index_dir, ref, gtf, threads, STAR):
    STARindexfile = ['chrLength.txt', 'chrName.txt', 'chrNameLength.txt', 'chrStart.txt', 'exonGeTrInfo.tab',
                     'exonInfo.tab', 'geneInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt',
                     'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']
    for STARfile in STARindexfile:
        if not os.path.isfile(os.path.join(star_index_dir, STARfile)):
            STAR_index_cmd = '%s --runMode genomeGenerate --genomeDir %s --runThreadN %s ' \
                             '--genomeFastaFiles %s --sjdbGTFfile %s' % (STAR, star_index_dir, threads, ref, gtf)
            if 'STAR_index' in myconf:
                STAR_index_cmd += join_parameter('STAR_index', myconf['STAR_index'])
            shell_script += STAR_index_cmd + '\n'
            break
    return shell_script

## fastq2bam，测序数据映射（bwa2）
def fastq2bam(myconf, bwa2, samtools, threads, ref, samindex, name, fastq, outdir, shell_script):
    map_cmd = '%s mem -t %s ' % (bwa2, threads)
    sammark = '-R "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s"' % (samindex, name)
    map_cmd += sammark + ' '
    if 'map' in myconf:
        map_cmd += join_parameter('bwa-mem2 mem', myconf['map'])
    sam, bam = os.path.join(outdir, name + '.sam'), os.path.join(outdir, name + '.bam')
    map_cmd += ref + ' ' + fastq[0] + ' ' + fastq[1] + ' -o %s\n' % (sam)
    map_cmd += '%s view -t %s -b -O bam -o %s %s\n' % (samtools, threads, bam, sam)
    map_cmd += 'rm -rf %s\n' % (sam)
    shell_script += map_cmd
    return shell_script, bam

## fastq2bam，测序数据映射（转录组，STAR）
def fastq2bam_STAR(myconf, STAR, samtools, star_index_dir, outdir, fastq, name, threads, shell_script):
    abs_prefix = os.path.join(outdir, name+'.')
    STAR_cmd = '%s --runThreadN %s --genomeDir %s --readFilesIn %s %s --outFileNamePrefix %s ' \
               '--outSAMtype BAM Unsorted' % (STAR, threads, star_index_dir, fastq[0], fastq[1], abs_prefix)
    if 'map_STAR' in myconf:
        STAR_cmd += join_parameter('map_STAR', myconf['map_STAR'])
    shell_script += STAR_cmd + '\n'
    star_bam = abs_prefix + 'Aligned.out.bam'
    shell_script, star_sbam = sort(myconf, samtools, threads, star_bam, outdir, shell_script)
    return shell_script, star_sbam

## BAM文件排序
def sort(myconf, samtools, threads, bam, outdir, shell_script):
    abs_prefix = os.path.basename(os.path.abspath(bam)).strip('.bam')
    sbam = os.path.join(outdir, abs_prefix + '.sorted.bam')
    sort_cmd = '%s sort -@ %s ' % (samtools, threads)
    if 'sort' in myconf:
        sort_cmd += join_parameter('samtools sort', myconf['sort'])
    sort_cmd += '%s > %s\n' % (bam, sbam)
    shell_script += sort_cmd
    return shell_script, sbam

## BAM文件去PCR重复
def markdup(myconf, sambamba, samtools, threads, bam, outdir, shell_script):
    abs_prefix = os.path.basename(os.path.abspath(bam)).strip('.bam')
    mbam = os.path.join(outdir, abs_prefix + '.markdup.bam')
    markdup_cmd = '%s markdup -t %s ' % (sambamba, threads)
    if 'markdup' in myconf:
        markdup_cmd += join_parameter('sambamba markdup', myconf['markdup'])
    markdup_cmd += '%s %s\n' % (bam, mbam)
    markdup_cmd += '%s index -@ %s %s\n' % (samtools, threads, mbam)
    shell_script += markdup_cmd
    return shell_script, mbam

## BAM文件还原为fastq
def bam2fastq(myconf, samtools, fastp, threads, bam, outdir, shell_script):
    abs_prefix = os.path.basename(os.path.abspath(bam)).strip('.bam')
    split_fastq = [os.path.join(outdir,abs_prefix+'_1.split.fastq'),os.path.join(outdir,abs_prefix+'_2.split.fastq')]
    split_cmd = '%s fastq -1 %s -2 %s -@ %s ' % (samtools, split_fastq[0], split_fastq[1], threads)
    if 'bam2fastq' in myconf:
        split_cmd += join_parameter('samtools fastq', myconf['bam2fastq'])
    split_cmd += '%s\n' % (bam)
    shell_script += split_cmd
    fix_fastq = [os.path.join(outdir, abs_prefix+'_1.split.fix.fastq'),os.path.join(outdir, abs_prefix+'_2.split.fix.fastq')]
    fix_html, fix_json = os.path.join(outdir, 'fix.html'), os.path.join(outdir, 'fix.json')
    fix_cmd = '%s -I %s -O %s -i %s -o %s -h %s -j %s\n' % (fastp, split_fastq[0], fix_fastq[0], split_fastq[1], fix_fastq[1], fix_html, fix_json)
    fix_cmd += 'rm -rf %s %s %s %s\n' % (fix_json, fix_html, split_fastq[0], split_fastq[1])
    shell_script += fix_cmd
    return shell_script, fix_fastq

"""使用gatk进行的变异识别"""
def variantcalling(bam, myconf, gatk, outdir, name, ref, shell_script):
    vcf = os.path.join(outdir, name + '.vcf')
    calling_cmd = '%s --java-options "%s" HaplotypeCaller ' % (gatk, '-Xmx' + myconf['java-options']['memory'] + 'g')
    calling_cmd += '-R %s -I %s -O %s' % (ref, bam, vcf)
    if 'call' in myconf:
        calling_cmd += join_parameter('gatk HaplotypeCaller', myconf['call'])
    shell_script += calling_cmd + '\n'
    return shell_script, vcf

"""使用featureCounts计数reads"""
def featureCounts(bam, myconf, featureCounts, gtf, ref, outdir, name, shell_script):
    readcount = os.path.join(outdir, name + '.featurecounts.txt')
    fc_cmd = '%s -a %s -G %s -F GTF -o %s ' % (featureCounts, gtf, ref, readcount)
    if 'featureCounts' in myconf:
        fc_cmd += join_parameter('featureCounts', myconf['featureCounts'])
    fc_cmd += '%s\n' % (bam)
    shell_script += fc_cmd
    return shell_script, readcount