import os
import sys
import numpy as np
import pandas as pd
from numpy.linalg import eig
from sklearn import manifold
from sklearn import preprocessing
import matplotlib.pyplot as plt

def PCA(features,  k, outdir, name, label_list, header_list, tg_label=None, threshold=0.95, decentral=False, maxabs=False, cov=False, var_plt=False, pca2D_plt=False, pca3D_plt=False):
    ## 标准化
    features, k1 = features, k
    if decentral==True:
        features = features - features.mean(axis=0)
    if maxabs == True:
        features = preprocessing.maxabs_scale(features, axis=0,copy=True)
    ## 计算协方差矩阵
    features_cov = np.cov(features.T, ddof=0)
    ## 计算特征值和特征向量
    evalues, evectors = eig(features_cov)
    evalues, evectors = np.real(evalues), np.real(evectors)
    ## 计算特征向量（主成分）的贡献率和累计贡献率
    tot = sum(evalues)
    var_exp = [(i/tot) for i in sorted(evalues, reverse=True)]
    cum_var_exp = np.cumsum(var_exp)
    ## 根据累计贡献率阈值确定取用的主成分个数
    k2 = list(cum_var_exp).index([x for x in cum_var_exp if x>=threshold][0])+1
    ku = max(k1, k2)
    ## 选取特征值最大的ku个特征向量
    klarge_index = evalues.argsort()[-ku:][::-1]
    k_evectors = evectors[klarge_index]
    ## 计算降维后的数据
    af_ddr = np.dot(features, k_evectors.T)
    ## 贡献率累计贡献率阈作图
    if var_plt==True:
        plt.bar(range(1,ku+1), var_exp[:ku], alpha=0.5, align='center',label='individual var')
        plt.step(range(1,ku+1), cum_var_exp[:ku], where='mid', label='cumulative var')
        plt.ylabel('Variance Ratio')
        plt.xlabel('Principal Components')
        plt.legend(loc='best')
        plt.savefig(os.path.join(outdir, name + '.PCA.variance.png'))
        plt.close()
    ## 将降维后的数据和特征向量保存为txt文件
    ku_order = ['PC'+str(x) for x in list(range(1, ku+1))]
    array2txt(af_ddr, os.path.join(outdir, name+'.PCA.ddr.txt'), header_list=ku_order, label_list=label_list)
    array2txt(k_evectors, os.path.join(outdir, name+'.PCA.vec.txt'), header_list=header_list, label_list=ku_order)
    exp_order = ['PC'+str(x) for x in list(range(1, len(header_list)+1))]
    array2txt(np.array(list(zip(var_exp, cum_var_exp))), os.path.join(outdir, name+'.PCA.exp.txt'), header_list=['Variance', 'Cumulative Variance'], label_list=exp_order)
    if cov == True:
        array2txt(features_cov, os.path.join(outdir, name+'.PCA.cov.txt'), header_list=label_list, label_list=label_list)
    with open(os.path.join(outdir,name+'.PCA.summary.txt'), 'w') as output:
        output.write('去中心化 = %s\n比例缩放 = %s\n' % (str(decentral), str(maxabs)))
        output.write('贡献率阈值 = %s\n预估主成分数 = %s\n实际主成分数 = %s\n' % (threshold, k1, k2))
        output.write('文件说明:\n*.PCA.ddr.txt = 降维后的数据\n')
        output.write('*.PCA.vec.txt = 主成分的向量拆解\n')
        output.write('*.PCA.exp.txt = 主成分贡献率及累计贡献率\n')
        if cov==True:
            output.write('*.PCA.cov.txt = 协方差矩阵\n')
    ## PCA作图
    if pca2D_plt==True and ku>= 2:
        Axes2Dplt(af_ddr, tg_label, outdir, name, 'PCA')
    if pca3D_plt==True and ku >= 3:
        Axes3Dplt(af_ddr, tg_label, outdir, name, 'PCA')

def array2txt(array, txt_path, header_list, label_list):
    df = pd.DataFrame(array, index=label_list, columns=header_list)
    df.to_csv(txt_path, sep='\t')

def Axes2Dplt(af_ddr, label_list, outdir, name, method):
    x = list(af_ddr[:,0]); y = list(af_ddr[:,1])
    plt.scatter(x,y,color='pink')
    if not label_list == None:
        for px,py,pl in zip(x,y,label_list):
            plt.text(px,py,pl,verticalalignment='center',horizontalalignment='right')
    plt.xlabel('%s 1' % (method.strip('A'))); plt.ylabel('%s 2' % (method.strip('A')))
    plt.savefig(os.path.join(outdir, name + '.'+ method + '.2D.png'))
    plt.close()

def Axes3Dplt(af_ddr, label_list, outdir, name, method):
    x = list(af_ddr[:,0]); y = list(af_ddr[:,1]); z = list(af_ddr[:,2])
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(x,y,z,color='red')
    if not label_list == None:
        for px,py,pz,pl in zip(x,y,z,label_list):
            ax.text(px,py,pz,pl,verticalalignment='center',horizontalalignment='right')
    ax.set_xlabel('%s 1' % (method.strip('A'))); ax.set_ylabel('%s 2' % (method.strip('A'))); ax.set_zlabel('%s 3' % (method.strip('A')))
    plt.savefig(os.path.join(outdir, name + '.'+ method + '.3D.png'))
    plt.close()

def t_SNE(features, outdir, name, label_list, tg_label=None, decentral=False, maxabs=False, ts2D=True, ts3D=False):
    ## 标准化
    if decentral == True:
        features = features - features.mean(axis=0)
    if maxabs == True:
        features = preprocessing.maxabs_scale(features, axis=0, copy=True)
    ## 降维2D
    if ts2D == True:
        ts_2D = manifold.TSNE(n_components=2, init='pca', random_state=0)
        features_ts_2D = ts_2D.fit_transform(features)
        Axes2Dplt(features_ts_2D, tg_label, outdir, name, 't-SNE')
        array2txt(features_ts_2D, os.path.join(outdir, name + '.t-SNE.2D.ddr.txt'), header_list=['t-SNE1', 't-SNE2'], label_list=label_list)
    ## 降维3D
    if ts3D == True:
        ts_3D = manifold.TSNE(n_components=3, init='pca', random_state=0)
        features_ts_3D = ts_3D.fit_transform(features)
        Axes3Dplt(features_ts_3D, tg_label, outdir, name, 't-SNE')
        array2txt(features_ts_3D, os.path.join(outdir, name + '.t-SNE.3D.ddr.txt'), header_list=['t-SNE1', 't-SNE2', 't-SNE3'], label_list=label_list)
    if ts2D == True or ts3D == True:
        with open(os.path.join(outdir,name+'.t-SNE.summary.txt'), 'w') as output:
            output.write('去中心化 = %s\n比例缩放 = %s\n' % (str(decentral), str(maxabs)))
            output.write('文件说明:\n')
            if ts2D == True:
                output.write('*.t-SNE.2D.ddr.txt = 2维降维后的数据\n')
            if ts3D == True:
                output.write('*.t-SNE.3D.ddr.txt = 3维降维后的数据\n')

def feature_read(feature_txt, label=True, tg_sample=None):
    feature_df = pd.read_csv(feature_txt, sep='\t', index_col=0, header=0)
    label_list = list(feature_df.index)
    tg_label =None
    if label == True:
        tg_label = list(feature_df.index)
        if tg_sample != None:
            tg_sample_list = tg_sample.split(',')
            tg_label = [None]*len(label_list)
            for tg in tg_sample_list:
                tg_label[int(tg)-1] = label_list[int(tg)-1]
    header_list = list(feature_df.columns)
    features = np.array(feature_df)
    return features, label_list, header_list, tg_label
