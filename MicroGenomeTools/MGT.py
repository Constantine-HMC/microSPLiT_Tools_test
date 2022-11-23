"""MGT(MicroGenomeTools)"""

"""
authors:
Zhuochong Liu, Flavius Constantine, Motochika Hata

version:
test.0.0.1

description:
    This is MicroGenomeTools, a python module for generating shell script of different common analysis
workflow with pre-set tools and pipeline, including variant caling (SNP and InDel), structural variant
detection, genome assembly and annotation, and so on. In addition, MicroGenomeTools provides a module 
to create your own workflow.
    Thanks for using MicroGenomeTools and wish you all the best.
"""

import sys
import os
import configparser
sys.path.append(sys.path[0])
import MGT_module

"""读取参数文件"""
def read_conf(conftxt):
    abs_path = os.path.abspath(conftxt)
    if not os.path.isfile(abs_path):
        raise FileNotFoundError('No such file: config: %s' % (abs_path))
    else:
        myconf = configparser.ConfigParser(strict=False)
        myconf.optionxform = str
        myconf.read(abs_path)
        return myconf

"""VariantCalling (SNP and InDel)"""
def variantcalling(inputtype, inputfile, ref, myconf, name, samindex, outdir, threads='4', cleanDir=False,
                   fp_TF=False, fqc_TF=False, sort_tf=False, markdup_tf=False):
    ## 空白脚本
    shell_script = '#!/usr/bin/bash\n'
    ## check file检查输入数据
    ### sequencing data测序数据文件
    input_list = inputfile.strip().split(',')
    for file in input_list:
        MGT_module.check_file(file)
    ### 参考基因组
    MGT_module.check_file(ref)
    ## check software检查所需软件
    fastqc, fastp = myconf['software']['fastqc'], myconf['software']['fastp']
    bwa2, samtools =  myconf['software']['bwa2'], myconf['software']['samtools']
    sambamba, gatk = myconf['software']['sambamba'], myconf['software']['gatk']
    MGT_module.check_software([fastqc, fastp, bwa2, samtools, sambamba, gatk])
    ## check index检查参考基因组的索引文件
    shell_script = MGT_module.ref_index(shell_script, ref, bwa2, samtools)
    ## 检查参考基因组的系列字典文件是否存在
    shell_script = MGT_module.ref_dict(shell_script, ref, gatk)
    ## check dir检查输出文件夹
    MGT_module.check_dir(outdir, True, cleanDir)
    ## 检查输入文件的类型
    if not inputtype in ['fastq', 'bam']:
        raise TypeError('Incorrect inputtype: %s, must be pair-end fastq or BAM' % (inputtype))
    ## 存储中间fastq、bam文件的文件夹
    fq_dir = os.path.join(outdir, 'fastq')
    if inputtype == 'fastq' and fp_TF == True:
        shell_script += 'mkdir %s\n' % (fq_dir)
    bam_dir = os.path.join(outdir, 'bam')
    bwa_bam_dir = os.path.join(bam_dir, 'bwa_bam')
    if (sort_tf == True or markdup_tf == True):
        shell_script += 'mkdir %s\n' % (bam_dir)
        shell_script += 'mkdir %s\n' % (bwa_bam_dir)
    ## fastqc质控及fastq2bam测序数据映射（bwa2）
    if inputtype == 'fastq':
        fastq = inputfile.strip().split(',')
        if fqc_TF == True:
            qc_result_dir = os.path.join(outdir, 'fastqc')
            shell_script += 'mkdir %s\n' % (qc_result_dir)
            row_dir = os.path.join(qc_result_dir, 'row')
            shell_script += 'mkdir %s\n' % (row_dir)
            shell_script = MGT_module.fq_fastqc(myconf, fastqc, fastq, row_dir, threads, shell_script)
        if fp_TF == True:
            qc_fq_dir = os.path.join(fq_dir, 'fastp')
            shell_script += 'mkdir %s\n' % (qc_fq_dir)
            shell_script, fastq = MGT_module.fq_fastp(myconf, fastp, name, fastq, qc_fq_dir, threads, shell_script)
            if fqc_TF == True:
                after_qc_dir = os.path.join(outdir, 'fastqc', 'fastp')
                shell_script += 'mkdir %s\n' % (after_qc_dir)
                shell_script = MGT_module.fq_fastqc(myconf, fastqc, fastq, after_qc_dir, threads, shell_script)
        shell_script, bam = MGT_module.fastq2bam(myconf, bwa2, samtools, threads, ref, samindex, name, fastq, bwa_bam_dir, shell_script)
    else:
        bam = inputfile
    ## sort
    sbam = bam
    if inputtype=='fastq' or (inputtype=='bam' and sort_tf==True):
        shell_script, sbam = MGT_module.sort(myconf, samtools, threads, bam, bwa_bam_dir, shell_script)
    ## markdup
    smbam = sbam
    if inputtype=='fastq' or (inputtype=='bam' and markdup_tf==True):
        shell_script, smbam = MGT_module.markdup(myconf, sambamba, samtools, threads, sbam, bwa_bam_dir, shell_script)
    ## variantcalling
    vcf_dir = os.path.join(outdir, 'vcf')
    shell_script += 'mkdir %s\n' % (vcf_dir)
    shell_script, vcf = MGT_module.variantcalling(smbam, myconf, gatk, vcf_dir, name, ref, shell_script)
    ## 写入脚本
    script = os.path.join(outdir, name+'.VariantCalling.sh')
    with open(script, 'w') as script_out:
        script_out.write(shell_script)
    return shell_script, bam, vcf

"""FeatureCounts"""
def FeatureCounts(inputtype, inputfile, ref, myconf, star_index_dir, gtf, name, samindex, outdir, threads='4',
                  cleanDir=False, fp_TF=False, fqc_TF=False, FQtype='untreated', BAMtype='bwa',
                  sort_tf=False, markdup_tf=False):
    ## 空白脚本
    shell_script = '#!/usr/bin/bash\n'
    ## check file检查输入数据
    ### sequencing data测序数据文件
    input_list = inputfile.strip().split(',')
    for file in input_list:
        MGT_module.check_file(file)
    ### 参考基因组及注释文件
    MGT_module.check_file(ref)
    MGT_module.check_file(gtf)
    ## check software检查所需软件
    fastqc, fastp = myconf['software']['fastqc'], myconf['software']['fastp']
    bwa2, samtools = myconf['software']['bwa2'], myconf['software']['samtools']
    sambamba = myconf['software']['sambamba']
    STAR, featureCounts = myconf['software']['STAR'], myconf['software']['featureCounts']
    MGT_module.check_software([fastqc, fastp, bwa2, samtools, sambamba, STAR, featureCounts])
    ## check index检查参考基因组的索引文件（bwa2）
    shell_script = MGT_module.ref_index(shell_script, ref, bwa2, samtools)
    ## check index检查参考基因组的索引文件（STAR）
    shell_script = MGT_module.ref_index_STAR(shell_script, myconf, star_index_dir, ref, gtf, threads, STAR)
    ## check dir检查输出文件夹
    MGT_module.check_dir(outdir, True, cleanDir)
    ## 检查输入文件类型
    if not inputtype in ['fastq', 'bam']:
        raise TypeError('Incorrect inputtype: %s, must be "fastq" or "BAM"' % (inputtype))
    if inputtype == 'bam':
        if not BAMtype in ['bwa','star']:
            raise TypeError('Incorrect input BAM file type: %s, must be "bwa" or "star"' % (inputtype))
    if inputtype == 'fastq':
        if not FQtype in ['untreated','treated']:
            raise TypeError('Incorrect input fastq file type: %s, must be "treated" or "untreated"' % (inputtype))
    ## 存储中间fastq、bam文件的文件夹
    fq_dir = os.path.join(outdir, 'fastq')
    if not (inputtype == 'bam' and BAMtype == 'star') and not (FQtype == 'treated'):
        shell_script += 'mkdir %s\n' % (fq_dir)
    bam_dir = os.path.join(outdir, 'bam')
    bwa_bam_dir = os.path.join(bam_dir, 'bwa_bam')
    if not (inputtype == 'bam' and BAMtype == 'star' and sort_tf == False):
        shell_script += 'mkdir %s\n' % (bam_dir)
        if not BAMtype == 'star' and not (sort_tf == False or markdup_tf == False):
            shell_script += 'mkdir %s\n' % (bwa_bam_dir)
    ## fastqc质控及fastq2bam测序数据映射（bwa2）
    fastq = inputfile.strip().split(','); bam = inputfile
    if inputtype == 'fastq' and FQtype == 'untreated':
        if fqc_TF == True:
            qc_result_dir = os.path.join(outdir, 'fastqc')
            shell_script += 'mkdir %s\n' % (qc_result_dir)
            row_dir = os.path.join(qc_result_dir, 'row')
            shell_script += 'mkdir %s\n' % (row_dir)
            shell_script = MGT_module.fq_fastqc(myconf, fastqc, fastq, row_dir, threads, shell_script)
        if fp_TF == True:
            qc_fq_dir = os.path.join(fq_dir, 'fastp')
            shell_script += 'mkdir %s\n' % (qc_fq_dir)
            shell_script, fastq = MGT_module.fq_fastp(myconf, fastp, name, fastq, qc_fq_dir, threads, shell_script)
            if fqc_TF == True:
                after_qc_dir = os.path.join(outdir, 'fastqc', 'fastp')
                shell_script += 'mkdir %s\n' % (after_qc_dir)
                shell_script = MGT_module.fq_fastqc(myconf, fastqc, fastq, after_qc_dir, threads, shell_script)
        shell_script, bam = MGT_module.fastq2bam(myconf, bwa2, samtools, threads, ref, samindex, name, fastq, bwa_bam_dir, shell_script)
    ## sort
    sbam = bam
    if (inputtype == 'fastq' and FQtype == 'untreated') or (inputtype == 'bam' and BAMtype == 'bwa' and sort_tf == True):
        shell_script, sbam = MGT_module.sort(myconf, samtools, threads, bam, bwa_bam_dir, shell_script)
    ## markdup
    smbam = sbam
    if (inputtype == 'fastq' and FQtype == 'untreated') or (inputtype == 'bam' and BAMtype == 'bwa' and markdup_tf == True):
        shell_script, smbam = MGT_module.markdup(myconf, sambamba, samtools, threads, sbam, bwa_bam_dir, shell_script)
    ## BAM文件还原为fastq
    if (inputtype == 'fastq' and FQtype == 'untreated') or (inputtype == 'bam' and BAMtype == 'bwa'):
        split_fq_dir = os.path.join(fq_dir,'split')
        shell_script += 'mkdir %s\n' % (split_fq_dir)
        shell_script, fastq = MGT_module.bam2fastq(myconf, samtools, fastp, threads, smbam, split_fq_dir, shell_script)
        if fqc_TF == True:
            split_qc_dir = os.path.join(outdir, 'fastqc', 'split')
            shell_script += 'mkdir %s\n' % (split_qc_dir)
            shell_script = MGT_module.fq_fastqc(myconf, fastqc, fastq, split_qc_dir, threads, shell_script)
    ## STAR比对
    STAR_outdir = os.path.join(bam_dir, 'star_bam')
    if not (inputtype == 'bam' and BAMtype == 'star' and sort_tf == False):
        shell_script += 'mkdir %s\n' % (STAR_outdir)
    if inputtype == 'fastq' or (inputtype == 'bam' and BAMtype == 'bwa'):
        shell_script, star_sbam = MGT_module.fastq2bam_STAR(myconf, STAR, samtools, star_index_dir, STAR_outdir, fastq, name, threads, shell_script)
    else:
        star_sbam = bam
        if sort_tf == True:
            shell_script, star_sbam = MGT_module.sort(myconf, samtools, threads, star_sbam, STAR_outdir, shell_script)
    ## 使用featureCounts计数reads
    featureCounts_dir = os.path.join(outdir,'fectureCounts')
    shell_script += 'mkdir %s\n' % (featureCounts_dir)
    shell_script, readcount = MGT_module.featureCounts(star_sbam, myconf, featureCounts, gtf, ref, featureCounts_dir, name, shell_script)
    ## 写入脚本
    script = os.path.join(outdir, name + '.FeatureCounts.sh')
    with open(script, 'w') as script_out:
        script_out.write(shell_script)
    return shell_script, readcount