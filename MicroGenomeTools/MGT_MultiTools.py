"""MGT_MultiTools(MicroGenomeTools_MultiTools)"""

"""
authors:
Zhuochong Liu, Flavius Constantine, Motochika Hata

version:
test.0.0.1

description:
    This is MGT_MultiTools, function module of creating shell scripts, or performing joint analysis for
multi samples, such as joint calling, and differential expression analysis.
"""

import sys
import os
import configparser
sys.path.append(sys.path[0])
import MGT
import MGT_module

"""VariantCalling of multi samples多样本变异识别"""
def multi_variantcalling(sample_list, outdir, create=False, threads='4'):
    abs_out_path = os.path.abspath(outdir)
    MGT_module.check_file(sample_list)
    MGT_module.check_dir(abs_out_path, create)
    sam_results = {}
    input_list = open(sample_list, 'r')
    scripts_list = open(os.path.join(abs_out_path, 'shell.list'), 'w')
    for line in input_list:
        if not line.startswith('#'):
            ll = line.strip().split('\t')
            name, samindex, inputtype, inputfile, ref, conf = ll[0], ll[1], ll[2], ll[3], ll[4], ll[5]
            create, clean, fp_TF, fqc_TF = ll[6], ll[7], ll[8], ll[9]
            sam_outdir = os.path.join(outdir, name)
            myconf = MGT.read_conf(conf)
            shell_script, bam, vcf = MGT.variantcalling(inputtype, inputfile, ref, myconf, name, samindex, sam_outdir, threads, create, clean, fp_TF, fqc_TF)
            sam_results[name]={'script':shell_script, 'bam':bam, 'vcf':vcf}
            scripts_list.write(shell_script+'\n')
    input_list.close(); scripts_list.close()
    return sam_results