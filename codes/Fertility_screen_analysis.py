
import pandas as pd
import numpy as np
import pickle
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
plot_folder='/Users/vikash/Documents/Projects/Claire/Fertility_screen/Figures/'

def changeHeaderFromSequencer():

    ## we will import two header files and change
    fertility_df=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/Fertility_screen_manifest_pool1.txt",sep='\t')
    ngi_df=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/NGI_sample_id_to_sample_name.txt",sep='\t')

    ## read barcode count file

    count_df=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/result_240120_barcode_counts_table.txt",sep='\t')
    input=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/input_vector.txt', sep='\t')

    tmp_df=count_df.loc[:,count_df.columns[2:]].copy()

    tmp_df.loc['col_mean',:]=tmp_df.sum(axis=0)
    tmp_df.loc[:,'row_mean']=tmp_df.sum(axis=1)
    ## required columns
    req_cols=tmp_df.columns[tmp_df.loc['col_mean']>5].to_list()[:-1]  ####
    req_rows= tmp_df.index[tmp_df['row_mean']>20].to_list()[:-1]       ####

    cols=['Gene',  'Barcodes']

    for col in req_cols:
        cols.append(col)

    ########## zeros

    ##########

    req_df=count_df.loc[req_rows,cols].copy()


    req_df['pbanka_id']=[item[0] for item in req_df['Gene'].str.findall('PBANKA_[0-9]*')]

    long_genes=req_df['Gene'].to_list()

    dict1={}
    l2=[]
    for g in input['gene_id'].to_list():
        for t in long_genes:
            if g in t :
                dict1[t]=g


    for g in input['gene_id'].to_list():

        if g not in dict1.values():
            l2.append(g)
    ## this code is ngi to sample

    ngi_to_sample=dict(zip( ngi_df['NGI Sample ID'],ngi_df['Your sample name']))
    sample_to_sample_des=dict(zip( fertility_df['Sample name'],fertility_df['Sample description 2']))

    for ngi,sam in ngi_to_sample.items():
        if sam in sample_to_sample_des.keys():
           req_df.columns=req_df.columns.str.replace(ngi,sample_to_sample_des[sam])

    req_df.columns=req_df.columns.str.replace('_F$','')
    req_df.columns=req_df.columns.str.replace('_R$','')
    ### filter and make average


    df_modi=req_df.copy()
    df_modi.set_index('pbanka_id',inplace=True)
    df_modi=df_modi.drop(columns=['Gene',  'Barcodes'])

    # save Original columns
    save_df=df_modi.copy()
    saveColumns=save_df.columns

    # we are going to remove sequencing samples and replace by '_R1' and '_R2'
    newcols=[]
    for c in saveColumns:
        if '_R1' in c.split('S')[1]:
            newcols.append(c.split('S')[0]+'R1')
        elif '_R2' in c.split('S')[1]:
            newcols.append(c.split('S')[0]+'R2')
        else:
            continue

    save_df.columns=newcols
    # remove PbSTM168_
    # save_df.columns= save_df.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    saveColumns=save_df.columns

    df_modi=save_df.copy()
    # Modified columns such that latter on we can map whih two columns are pair end sequencing
    pair_read_columns = df_modi.columns.str.replace("_R[12]", '', regex=True)
    df_modi.columns =  pair_read_columns
    pair_read_map=getGroupMap(saveColumns, df_modi.columns)

    # Modified columns such that latter on we could analyze PCR error
    df_modi.columns= df_modi.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    pcr_columns = df_modi.columns.str.replace("_r[12]", '', regex=True)
    df_modi.columns = pcr_columns
    map_pcr=getGroupMap(saveColumns, df_modi.columns)

    # Modified columns such that latter on we could analyze disection error
    # to find mf1.1 vs mf1.2 and  mf2.1 vs mf2.2
    mossifeed_df=df_modi.copy()
    mf_df=df_modi.copy()
    # midgut_flag = df_modi.columns.str.extract("^[0-9]+_(.+mf[12].[12])$")
    midgut_flag = df_modi.columns.str.extract("_(.+mf[12].[12])$")
    midgut_columns = saveColumns[np.where(~midgut_flag.isnull())[0]]
    # dissection error
    midgut_df = save_df.loc[:, midgut_columns]
    midgut_df.columns= midgut_df.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    midgut_df.columns = midgut_df.columns.str.replace("_r[12]_R[12]$", '', regex=True)
    midgut_df.columns = midgut_df.columns.str.replace(".[12]$", '', regex=True)

    map_midgut = getGroupMap(midgut_columns, midgut_df.columns)

    # find mosquito feed diffrence mf1 and mf2 mosquito feed

    mossifeed_flag = mossifeed_df.columns.str.extract("_(.+mf[12])+")
    mossifeed_columns = saveColumns[np.where(~ mossifeed_flag.isnull())[0]]
    mossifeed_df = save_df.loc[:, mossifeed_columns]

    #

    mossifeed_df.columns= mossifeed_df.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    mossifeed_df.columns = mossifeed_df.columns.str.replace("_r[12]_R[12]$", '', regex=True)
    # colNames1=mossifeed_df.columns.str.replace("_R[12]$", '',regex=True)
    # mossifeed_df.columns=colNames1
    # colNames2=mossifeed_df.columns.str.replace("_r[12]", '',regex=True)
    # mossifeed_df.columns=colNames2
    colNames4=mossifeed_df.columns.str.replace("\\.[12]$", '', regex=True)
    mossifeed_df.columns=colNames4
    colNames5=mossifeed_df.columns.str.replace("\\_mf[12]$", '', regex=True)
    mossifeed_df.columns=colNames5

    map_mossifeed = getGroupMap(mossifeed_columns, mossifeed_df.columns)

    ### group  to plot mosquito feed 1 and mosquito feed 2

    map_mf1=find_mf1_columns(mf_df,saveColumns,save_df)
    map_mf2=find_mf2_columns(mf_df,saveColumns,save_df)


    req_df.to_csv('/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/Filtered_Count.txt',sep='\t',index=None)

    return req_df,save_df,map_pcr,map_midgut,map_mossifeed,map_mf1,map_mf2


def find_mf1_columns(mossifeed_df,saveColumns,save_df):

    # find mosquito feed diffrence mf1 and mf2 mosquito feed

    mossifeed_flag = mossifeed_df.columns.str.extract("_(.+mf[1])+")
    mossifeed_columns = saveColumns[np.where(~ mossifeed_flag.isnull())[0]]
    mossifeed_df = save_df.loc[:, mossifeed_columns]

    #

    mossifeed_df.columns= mossifeed_df.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    mossifeed_df.columns = mossifeed_df.columns.str.replace("_r[12]_R[12]$", '', regex=True)
    # colNames1=mossifeed_df.columns.str.replace("_R[12]$", '',regex=True)
    # mossifeed_df.columns=colNames1
    # colNames2=mossifeed_df.columns.str.replace("_r[12]", '',regex=True)
    # mossifeed_df.columns=colNames2
    colNames4=mossifeed_df.columns.str.replace("\\.[12]$", '', regex=True)
    mossifeed_df.columns=colNames4
    colNames5=mossifeed_df.columns.str.replace("\\_mf[1]$", '', regex=True)
    mossifeed_df.columns=colNames5

    map_mossifeed = getGroupMap(mossifeed_columns, mossifeed_df.columns)
    return map_mossifeed



def find_mf2_columns(mossifeed_df,saveColumns,save_df):

    # find mosquito feed diffrence mf1 and mf2 mosquito feed

    mossifeed_flag = mossifeed_df.columns.str.extract("_(.+mf[2])+")
    mossifeed_columns = saveColumns[np.where(~ mossifeed_flag.isnull())[0]]
    mossifeed_df = save_df.loc[:, mossifeed_columns]

    #

    mossifeed_df.columns= mossifeed_df.columns.str.replace("PbFertility1_[0-9]+_", '', regex=True)
    mossifeed_df.columns = mossifeed_df.columns.str.replace("_r[12]_R[12]$", '', regex=True)
    # colNames1=mossifeed_df.columns.str.replace("_R[12]$", '',regex=True)
    # mossifeed_df.columns=colNames1
    # colNames2=mossifeed_df.columns.str.replace("_r[12]", '',regex=True)
    # mossifeed_df.columns=colNames2
    colNames4=mossifeed_df.columns.str.replace("\\.[12]$", '', regex=True)
    mossifeed_df.columns=colNames4
    colNames5=mossifeed_df.columns.str.replace("\\_mf[2]$", '', regex=True)
    mossifeed_df.columns=colNames5

    map_mossifeed = getGroupMap(mossifeed_columns, mossifeed_df.columns)
    return map_mossifeed



def readInputandCountFile(df):
    input=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/input_vector.txt', sep='\t')

    ### systamatic tests
    # 1) gneral view of contamination we need to find not common genes
    contamin_genes=set(df.index.to_list())-set(input['gene_id'].to_list())
    if len(contamin_genes)>0:

        conta_df=df.loc[contamin_genes,:].copy()
        print('Contaminate genes are')

    else:
        print('There is no contaminate genes ')
    ## identify drop out genes
    # find at vitro level
    tmp=df.copy()
    saveColumns=tmp.columns
    vitro_flag = tmp.columns.str.extract("_(t[12])+")
    vitro_columns = saveColumns[np.where(~ vitro_flag.isnull())[0]]

    drop_df1= df.loc[:,vitro_columns].copy()
    drop_df1.loc[:,'row_mean']=drop_df1.sum(axis=1)
    ## required columns
    dp_genes1= drop_df1.index[drop_df1['row_mean']<50].to_list()[:-1]       ####

    ####  identify drop outs at input level (d0)

    input_flag = tmp.columns.str.extract("_(d[0])+")
    input_columns = saveColumns[np.where(~ input_flag.isnull())[0]]

    drop_df2= df.loc[:,input_columns].copy()
    drop_df2.loc[:,'row_mean']=drop_df2.sum(axis=1)
    ## required columns
    dp_genes2= drop_df2.index[drop_df2['row_mean']<50].to_list()[:-1]       ####

    ### real input genes which are not drop out

    real_input_genes=set(input['gene_id'].to_list())-(set(dp_genes2)|set(dp_genes1))

    ### drop out real
    drop_out=real_input_genes-set(df.index.to_list())
    real_input_genes=real_input_genes-drop_out

    req_df=df.loc[real_input_genes,:].copy() ### This is the require dataframes

    return req_df


def calRel(df):
    n=df.shape[1]
    df=df.div(df.sum(axis=1), axis=0)
    # get log value
    # # cut=1e-10 # this is the cutoff
    # df[df<cut]=cut
    df=np.log2(df)
    #df=np.log10(df)
    return df



def cal_mean_and_sd_groupby_columns(df_modi,mapping):
    tmp_mean_T = pd.DataFrame(index=df_modi.index)
    tmp_std_T = pd.DataFrame(index=df_modi.index)
    for k,item in mapping.items():

        tmp_mean_T[k]=df_modi[item].mean(axis=1).copy()
        tmp_std_T[k]=df_modi[item].std(axis=1).copy()
    return tmp_mean_T,tmp_std_T


def getColumnsFormDF(df,tag):
    t_cols=[]
    for t in tag:
        flag = df.columns.str.contains(t)
        cols= df.columns[np.where(flag)[0]]
        t_cols.append(cols)
    return t_cols

def calDissectionRatio(df_m,df_sd,time,time2,mf,backgrounds):
    # based on pcr mean and sd we will calculate ratios
    # df_sd = df_sd / df_m  # this is just to take only mutants
    # time= ['d0', 'd7', 'd14']
    # mf=['mf1','mf2']
    # backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']


    ratioDict={}
    for b in backgrounds:
        for t in time:
            for f in mf:
                string=b+'_'+t+'_'+f
                time_cols = getColumnsFormDF(df_m, [string])

                ratioDict[(b,t,f)]=time_cols[0]

    time = time2
    rel_fit=pd.DataFrame(index=df_m.index)
    # rel_fit1 = pd.DataFrame(index=df_m.index)
    rel_err = pd.DataFrame(index=df_m.index)
    # rel_err1 = pd.DataFrame(index=df_m.index) #without log scale
    for b in backgrounds:
        for t in time:
            for f in mf:
                tcol = ratioDict[(b, 'd0', f)]
                input=pd.DataFrame(index=df_m.index)
                input_err = pd.DataFrame(index=df_sd.index)
                input = pd.concat([input, df_m[tcol].copy()], axis=1)
                input_err = pd.concat([input_err, df_sd[tcol].copy()], axis=1)  # error

                tcol=ratioDict[(b, t, f)]
                tech_df = pd.DataFrame(index=df_m.index)
                tech_df_err = pd.DataFrame(index=df_sd.index)
                tech_df = pd.concat([tech_df, df_m[tcol].copy()], axis=1)
                tech_df_err = pd.concat([tech_df_err, df_sd[tcol].copy()], axis=1)  # error


                rel_fit = pd.concat([rel_fit, tech_df.iloc[:, 0:].sub(input.iloc[:, 0], axis=0)], axis=1)  # error

                # SD on log2 scale
                # tech_df_err=np.log2(tech_df_err)
                # input_err=np.log2(input_err)

                tmp=tech_df_err.pow(2, axis=1).iloc[:,0:].add(input_err.pow(2,axis=1).iloc[:,0],axis=0)
                # rel_err=pd.concat([rel_err,tmp.apply(np.sqrt)],axis=1) SD
                rel_err=pd.concat([rel_err,tmp],axis=1) # varaiance




    return rel_fit, rel_err

def gaussianMeanAndVariance(rel_fit,rel_err):

    weight = rel_err.rdiv(1)

    nume_df=pd.DataFrame(rel_fit.values * weight.values, columns=rel_fit.columns, index=rel_fit.index)
    nume= nume_df.sum(axis=1)
    deno=weight.sum(axis=1)
    # calculate mean
    mean_df=nume/deno

    # calculate first varaiance
    var1= weight.sum(axis=1).rdiv(1)

    # claculate second variance
    term1=rel_fit.iloc[:, 0:].sub(mean_df, axis=0).pow(2,axis=1)
    term3=term1.div(rel_err,axis=0)
    term4=term3.sum(axis=1)
    term2=(1/(rel_fit.shape[1]-1))


    var2=var1*term4*term2

    var_max=pd.concat([var1, var2], axis=1).max(axis=1)
    sd_max=var_max.apply(np.sqrt)

    return [mean_df,sd_max,var_max]


def plotfitScatter(x,y,xlab,ylab,pdf,title):
    fig= plt.figure(figsize=(10,8))
    ax = plt.gca()
    xx=[]
    yy=[]
    pp=[]
    xpxp=[]


    for i in range(len(x)):
        xx.append(np.array(x[i]))
        yy.append(np.array(y[i]))
        z=np.polyfit(x[i], y[i], 1)
        p=np.poly1d(z)
        mini=np.min(x[i])
        maxi=np.max(x[i])
        pp.append(p)
        xp = np.linspace(mini, maxi, 100)
#         print(xp,pp)
        xpxp.append(xp)

        # # exponnetial fit
        # popt, pcov = curve_fit(exponenial_func, xx[i], yy[i], p0=(1, 1e-3, 1))
        # print(popt)
        # y3.append(exponenial_func(xp, *popt))
    lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], pp[0](xpxp[0]), 'r-',xx[1], yy[1], 'b.', xpxp[1], pp[1](xpxp[1]), 'b-',linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=15)
    plt.ylabel(ylab,fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(lines,['d0','d0','d13','d13'],fontsize=10);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);
#     plt.xlim([-14,1])
    plt.savefig(pdf)

    plt.show()


def calRel_withoutLog(df):

    n=df.shape[1]
    df=df.div(df.sum(axis=1), axis=0)
    # get log value
    # cut=1e-10
    # df[df<cut]=cut;
    return df


def plotfitScatterDay(x,y,xlab,ylab,pdf,title):
    fig= plt.figure(figsize=(10,8))
    ax = plt.gca()
    xx=[]
    yy=[]
    pp=[]
    xpxp=[]
    y3=[]
    for i in range(len(x)):
        xx.append(np.array(x[i]))
        yy.append(np.array(y[i]))
        z=np.polyfit(x[i], y[i], 1)
        p=np.poly1d(z)
        mini=np.min(x[i])
        maxi=np.max(x[i])
        pp.append(p)
        xp = np.linspace(mini, maxi, 100)
#         print(xp,pp)
        xpxp.append(xp)

        # # exponnetial fit
        # popt, pcov = curve_fit(exponenial_func, xx[i], yy[i], p0=(1, 1e-3, 1))
        # print(popt)
        # y3.append(exponenial_func(xp, *popt))
    lines=ax.plot(xx[0], yy[0], 'b.', xpxp[0], pp[0](xpxp[0]),linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=15)
    plt.ylabel(ylab,fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(lines,['d13','d13'],fontsize=10);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);

    plt.savefig(pdf)

    plt.show()



def calPvalZscore(cmb_fitness,geneMap):
    import scipy.stats as st

    import scipy.special as sp
    # anot, geneMap=readAnotation()
    # create df for ploting
    for k in cmb_fitness.keys():
        key1=k
        break
    comm=set(cmb_fitness[key1][0].index)

    geneName=[geneMap[item] for item in comm]
    viz_df=pd.DataFrame(index=cmb_fitness[key1][0].index)
    viz_df['geneName'] = geneName

    pval=0.01
    pval1=0.01
    upcut=np.log2(1)
    lowcut=np.log2(0.03)

    pheno_pval={}

    for key,item in cmb_fitness.items():

    # we will calculate p-vlaue and z score
        m=item[0]
        s=item[1]
        z=(upcut-m)/s
        df = pd.DataFrame(m, columns=[key])
        df[(key,'SD')]=s
        df['z_score']=z
        # df['pvalue']=1-sp.ndtr(df['z_score'])
        # df['pvalue_2tail'] = 2*(1 - sp.ndtr(df['z_score']))
        df['pvalue']=  (1 - st.norm.cdf(z))

       # we will calculate p value form 0.1
        z = (lowcut - m) / s
        df['z_score_2'] = z
        df['pvalue_2'] = (1 - st.norm.cdf(z))

        # select fast growing, NE, E, middle
        # if pvalue<0.05 and relative fitness >0
        df_fast=df[(df['pvalue'] < pval) & (df[key] > upcut)]
        df_ambiguous1=df[(df['pvalue'] < pval) & (df[key] < upcut)]
        NE_df=df[(df['pvalue'] > pval)]
        df_subset1= df[(df['pvalue_2'] >pval1)]
        df_ambiguous2 = df[(df['pvalue_2'] < pval1) & (df[key] > lowcut)]
        E_df=df[(df['pvalue_2'] < pval1)]

        # non essential genes
        NE=(NE_df.index)
        # essential genes
        E=E_df.index
        dropIdx=NE.union(E)
        # other genes
        df_other=df.drop(index=dropIdx)

        other=df_other.index

        viz_df[key[0]+'_'+key[1]+'_pheno'] = 'NA'
        viz_df[key[0]+'_'+key[1]+'_pheno'][E] = 'E'
        viz_df[key[0]+'_'+key[1]+'_pheno'][NE] = 'NE'
        viz_df[key[0]+'_'+key[1]+'_rel'] = m
        viz_df[key[0] + '_' + key[1] + '_sd'] = s

        pheno_pval['E']=E
        pheno_pval['NE'] = NE
        pheno_pval['other'] = other
        #import pdb;pdb.set_trace()  # we are checking the pvalues for normal genes .

    return viz_df



def plotMaleFemaleScatter(df):
    # fig = plt.figure(figsize=(8,8))
    fig, ax = plt.subplots(figsize=(9, 9))

    ###
    # change df_values

    df['145480_d13_pheno']=df['145480_d13_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
    df['GCKO2_d13_pheno']=df['GCKO2_d13_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
    # viz_df['Published_cross_phenotype']=viz_df['Published_cross_phenotype'].replace({'N': 'NA'})

    # cmap={'FM':"#66c2a5", 'RM':"#8da0cb", 'IM':"#fc8d62"}
    cmap_male={'FM':"#1b9e77", 'RM':"#7570b3", 'IM':"#d95f02"}
    cmap_female={'FF':"#1b9e77", 'RF':"#7570b3", 'IF':"#d95f02"}


    for i,item in enumerate(df['GCKO2_d13_pheno'].to_list()):
        marker_style = dict(color=cmap_female[item],markersize=6, markerfacecoloralt=cmap_male[df['145480_d13_pheno'][i]],markeredgecolor='white',alpha=0.8)
        ax.plot(df['GCKO2_d13_rel'][i],df['145480_d13_rel'][i], 'o',fillstyle='left', **marker_style)


    ax.set_xlabel('Relative growth rate (Female)',fontsize=15)
    ax.set_ylabel('Relative growth rate (Male)',fontsize=15)
    ax.tick_params(axis='y',labelsize=12)
    ax.tick_params(axis='x',labelsize=12)
    ax.set_xlim(-11,3)
    ax.set_ylim(-11,3)
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles, labels)
    infertile= mpatches.Patch(color="#fc8d62", label='Infertility (F/M)')
    fertile=mpatches.Patch(color="#66c2a5", label='Normal fertility (F/M)')
    reduced_fertile=mpatches.Patch(color="#8da0cb", label='Reduced fertility (F/M)')
    plt.legend(handles=[infertile,fertile,reduced_fertile],loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3, borderaxespad=0, frameon=False, prop={'size': 6})

    fig.savefig(plot_folder + "scatter_plot_male_female_RGR_pool1.pdf")




def plotInSnsNew(viz_df):



    # sns.set_style("white")
    sns.set(style="white",font_scale=1.5)
    print('hi')

    # change df_values

    viz_df['145480_d13_pheno']=viz_df['145480_d13_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
    viz_df['GCKO2_d13_pheno']=viz_df['GCKO2_d13_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
    # viz_df['Published_cross_phenotype']=viz_df['Published_cross_phenotype'].replace({'N': 'NA'})


    ##################### BAR PLOTS BEGINS


    plt.figure(figsize=(16*4, 20*4))
    bar_df=viz_df.sort_values(by=["145480_d13_rel"])

    # cmap={'FM':"#66c2a5", 'RM':"#8da0cb", 'IM':"#fc8d62"}
    cmap={'FM':"#1b9e77", 'RM':"#7570b3", 'IM':"#d95f02"}

    hue=[ "#1b9e77",  "#7570b3" ,"#d95f02"]
    clrs = [cmap[bar_df.loc[x,"145480_d13_pheno"]]  for x in bar_df.index]

    ax = sns.barplot(x="145480_d13_rel", y="geneName", data=bar_df, ci=None,palette = clrs)
    # new_labels = ['FM', 'RM' , 'IM']
    # for t, l in zip(ax._legend.texts, new_labels): t.set_text(l)

    plt.errorbar(x=bar_df["145480_d13_rel"], y=ax.get_yticks(), xerr=bar_df["145480_d13_sd"], fmt='none',c='gray')

    patch =[ mpatches.Patch(color="#1b9e77", label='FM'), mpatches.Patch(color="#7570b3", label='RM'), mpatches.Patch(color="#d95f02", label='IM')]
    plt.legend(handles=patch)
    fig = ax.get_figure()
    fig.savefig(plot_folder + "bar_145480_plot_RGR_pool1.pdf")



    plt.figure(figsize=(16*4, 20*4))
    bar_df=viz_df.sort_values(by=["GCKO2_d13_rel"])

    cmap={'FF':"#1b9e77", 'RF':"#7570b3", 'IF':"#d95f02"}


    clrs = [cmap[bar_df.loc[x,"GCKO2_d13_pheno"]]  for x in bar_df.index]

    ax = sns.barplot(x="GCKO2_d13_rel", y="geneName", data=bar_df, ci=None,palette = clrs)
    plt.errorbar(x=bar_df["GCKO2_d13_rel"], y=ax.get_yticks(), xerr=bar_df["GCKO2_d13_sd"], fmt='none',c='gray')
    patch =[ mpatches.Patch(color="#1b9e77", label='FF'), mpatches.Patch(color="#7570b3", label='RF'), mpatches.Patch(color="#d95f02", label='IF')]
    plt.legend(handles=patch)
    fig = ax.get_figure()
    fig.savefig(plot_folder + "bar_GCKO2_plot_RGR_pool1.pdf")


def stepByStep_barSeqAnalysis():

    ############# STEP1 ###########

    ### we are taking gene counts which we found from barseq experiment
    ## arranged columns based on diffrent kind of errors:
    ## midgut errors, disection errors, PCR errors
    df_req, df_data,map_pcr, map_midgut, map_mossifeed,map_mf1,map_mf2 = changeHeaderFromSequencer()


    ############# STEP 2 ###########

    # we want to find out where there is contamination because due to experiment it could be some mutanats can be mixed.
    # for this we need to read input file for pools and then we need to match with what we get from barseq data

    # read Input pool file and check with data

    df_ori=readInputandCountFile(df_data)


    ### ######  Step 3  calculate relative abundance
    df_ori = df_ori + 1 ### this trick we used to remove zero count
    ref_ori = calRel(df_ori.T) # calculate relative frequency
    ref_ori = ref_ori.T # tarnspose



    ##################### plot relative log2 counts ####################

    input_genes=ref_ori.index.to_list()
    ## map previous id
    geneConv=getNewIdfromPrevID(input_genes)

    ## mf1 mean and sd
    mf1_m, mf1_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mf1)
    mf2_m, mf2_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mf2)

    # ## plot relative fitnees
    # out_pdf = plot_folder+'relative_fitness.pdf'
    # plot_relativeFitness(mf1_m, mf1_sd,mf2_m, mf2_sd,geneConv,out_pdf)

    ### end for plotting relative fitness


    ######### we are going to plot errors

    #calculate error for PCR duplicates

    # df_ori=df_ori
    # import pdb;pdb.set_trace()
    ref_ori=calRel_withoutLog(df_ori.T)
    ref_ori=ref_ori.T
    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr)  # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed) # this is the mean and difference between mosiquito feeds

    tag=['d0_mf','d13_mf']

    ######## PCR error

    pcr_time_df,pcr_time_df_err,pcr_time_df_err1,pcr_time_df_nl,pcr_values,pcr_values_err,pcr_values_err1,pcr_values_nl=selectColumns(pcr_m, pcr_sd,tag)
    pdf=plot_folder+'PCR_error_relative abundance_pool1.pdf'

    #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    # plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','SD of relative abundances',pdf,'PCR_error')
    # plotfitScatter(pcr_values,pcr_values_err1,'log2(mean of relative abundances)','(SD of relative abundances)/mean',pdf,'PCR_error')

    ####
    tag=['d13_mf']
    gut_time_df,gut_time_df_err,gut_time_df_err1,gut_time_df_nl,gut_values,gut_values_err,gut_values_err1,gut_values_nl=selectColumns( midgut_m, midgut_sd,tag) # find only for mutants

    pdf=plot_folder+'dissection_error_pool1.pdf'
    #plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    # plotfitScatterDay(gut_values,gut_values_err1,'log2(mean relative abundances)','SD of relative abundances',pdf,'Dissection_error')


    ### we are going to plot distribution on dissection error plot

    # gut_time_df,gut_time_df_err1=plotHistogramonError(gut_values,gut_values_err1,pcr_time_df)



    #######   STEP4 calculate sctter plot male vs female


    ref_ori = calRel(df_ori.T) # calculate relative frequency
    ref_ori = ref_ori.T # tarnspose
    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr)  # this is the mean and difference between pcr duplicates

    # midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples
    #
    # mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed) # this is the mean and difference between mosiquito feeds



    #### calculate dissection error
    time= ['d0', 'd13']
    mf=['mf1','mf2']
    backgrounds=['GCKO2','145480']
    time2= ['d13']

    rel_fitness,rel_var= calDissectionRatio(pcr_m, pcr_sd,time,time2,mf,backgrounds)

    #

    # withlogtuple=tl[0]
    # nologtuple=tl[1]
    # normalized by reference poposed by oliver
    # control_genes=['PBANKA_102460','PBANKA_100210','PBANKA_100220','PBANKA_050120']

    control_genes=['PBANKA_102460' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']




    control_gene_info={}
    for col in rel_fitness.columns:
        control_fitness=rel_fitness.loc[control_genes,[col]].copy()
        control_var=rel_var.loc[control_genes,[col]].copy()
        l=gaussianMeanAndVariance(control_fitness.T,control_var.T)
        control_gene_info[col]=l # l[0] mean l[1]  SD   l[2] variance

    #  we want to normalize by control genes
    normalized_fit=rel_fitness.copy()
    normalized_var=rel_var.copy()

    for col in rel_fitness.columns:
        ctr_mean=control_gene_info[col][0][col]  #  0 mean
        ctr_var=control_gene_info[col][2][col]  # 1 variance

        # this is the relative mean
        normalized_fit.loc[:,col]=rel_fitness.loc[:,col].sub(ctr_mean)
        # relative variance on log scale

        normalized_var.loc[:,col]=rel_var.loc[:,col].add(ctr_var)


    cmb_fitness={}
    time= ['d13']
    backgrounds=['GCKO2','145480']

    for b in backgrounds:
        for t in time:
            string=b+'_'+t
            tcols = getColumnsFormDF(normalized_fit, [string])
            tm=normalized_fit[tcols[0]].copy()
            terr=normalized_var[tcols[0]].copy()

            cmb_fitness[(b,t)]=gaussianMeanAndVariance(tm,terr)

    # #plotrelFitness(cmb_fitness)
    # ##########  This function  is used for preparing plots for fertility website
    # res_df=prepareForPGfertility(cmb_fitness)


    viz_df=calPvalZscore(cmb_fitness,geneConv)



    ### create sctter plot male female relative fitness
    # plotMaleFemaleScatter(viz_df)



    plotInSnsNew(viz_df)

def getValues(df):
    df_values=[]
    for t_values in df.values:
        for val in t_values:
            df_values.append(val)
    return df_values

def getColumnsFormDF(df,tag):
    t_cols=[]
    for t in tag:
        flag = df.columns.str.contains(t)
        cols= df.columns[np.where(flag)[0]]
        t_cols.append(cols)
    return t_cols

def selectColumns(df_m, df_sd,tag):
    #with this function we want to compute error for each time point
    # df_m is the mean df_sd is SD.
    # test whether there is duplicates in genes
    if len(df_m.index) !=len(df_m.index.unique()):
        print ("there is duplicate genes: We should fixed which gene do we want to take in analysis")
        exit()


    time_cols=getColumnsFormDF(df_m,tag)
    df_rel=df_sd/df_m # this is just to take only mutants

    time_df_nolog=[]
    time_df=[]
    time_df_err=[]
    time_df_err1=[]
    # technical replicates

    for tcol in time_cols:
        tech_df = pd.DataFrame(index=df_m.index)
        tech_df_err = pd.DataFrame(index=df_sd.index)
        tech_df_err1 = pd.DataFrame(index=df_rel.index)
        tech_df = pd.concat([tech_df,df_m[tcol].copy()],axis=1)
        tech_df_err = pd.concat([tech_df_err,df_sd[tcol].copy()],axis=1) #  error
        tech_df_err1 = pd.concat([tech_df_err1,df_rel[tcol].copy()],axis=1) #  error
        cut=1e-10
        if tech_df.isnull().values.any():
            print ("relative abundance contains nan value do you want to ignore? we put 1")
            tech_df=tech_df.fillna(1)
        if tech_df_err.isnull().values.any():
            print ("relative abundance error contains nan value do you want to ignore? we put 0")
            tech_df_err=tech_df_err.fillna(0)
        # import pdb;pdb.set_trace()
        # tech_df[tech_df < cut]=cut
        time_df_nolog.append(tech_df)
        tech_df=np.log2(tech_df)
        # tech_df_err[tech_df_err < cut]=cut
        # tech_df_err=np.log2(tech_df_err)

        time_df.append(tech_df)
        time_df_err.append(tech_df_err)
        time_df_err1.append(tech_df_err1)

    # collaspe in one vector calculate errors and values
    l_values=[]
    for item in time_df:
# technical replicates ratio
        values=getValues(item)
        l_values.append(values)

    l_values_no_log = []
    for item in time_df_nolog:
        values = getValues(item)
        l_values_no_log.append(values)

    l_values_err=[]
    for item in time_df_err:
# technical replicates error
        values=getValues(item)
        l_values_err.append(values)

    l_values_err1=[]
    for item in time_df_err1:
# technical replicates error
        values=getValues(item)
        l_values_err1.append(values)

    return time_df,time_df_err,time_df_err1,time_df_nolog,l_values,l_values_err,l_values_err1,l_values_no_log

def plot_relativeFitness(mf1_m, mf1_sd,mf2_m, mf2_sd,geneConv,out_pdf):
    pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
    n=2
    genes=mf1_m.index;

    labels=[]
    for gene in genes:
        if gene in geneConv.keys():
            labels.append(geneConv[gene])
        else:
            labels.append(gene)

    ### get number of subplots

    num_subplots=2*len(genes) ### number of subplots

    l=len(genes)
    loops=int(l/(n*n))
    rem=l % (n*n)
    k=0

    per_page=num_subplots/16

    for loop in range(loops):
        fig = plt.figure(figsize=(15,15))
        for i in range(1,(n*n)+1):

            plt.subplot(n, n, i)
            x=[1,2]


            ### mf1 female
            y=[mf1_m.loc[genes[k],'GCKO2_d0'],mf1_m.loc[genes[k],'GCKO2_d13']]
            yerr=[mf1_sd.loc[genes[k],'GCKO2_d0'],mf1_sd.loc[genes[k],'GCKO2_d13']]
            plt.errorbar(x,y, yerr=yerr, fmt='b--')

            ### mf2 female
            y=[mf2_m.loc[genes[k],'GCKO2_d0'],mf2_m.loc[genes[k],'GCKO2_d13']]
            yerr=[mf2_sd.loc[genes[k],'GCKO2_d0'],mf2_sd.loc[genes[k],'GCKO2_d13']]
            plt.errorbar(x,y, yerr=yerr, fmt='r--')

            ### mf1 male
            y=[mf1_m.loc[genes[k],'145480_d0'],mf1_m.loc[genes[k],'145480_d13']]
            yerr=[mf1_sd.loc[genes[k],'145480_d0'],mf1_sd.loc[genes[k],'145480_d13']]
            plt.errorbar(x,y, yerr=yerr, fmt='b-')

            ### mf2 male
            y=[mf2_m.loc[genes[k],'145480_d0'],mf2_m.loc[genes[k],'145480_d13']]
            yerr=[mf2_sd.loc[genes[k],'145480_d0'],mf2_sd.loc[genes[k],'145480_d13']]
            plt.errorbar(x,y, yerr=yerr, fmt='r-')

            plt.ylabel('log2 relative fitness')
            plt.title(labels[k])
            plt.legend(('mf1_female', 'mf2_female','mf1_male', 'mf2_male'))
            plt.xticks([1, 1.5, 2],['day0', '', 'day13'],fontsize=15)

            plt.ylim(-18, 1)

            k=k+1
        pdf.savefig(fig)

    ## for the remaing one
    fig = plt.figure(figsize=(15,15))
    for i in range(1,rem+1):
        fig = plt.figure(figsize=(15,15))


        plt.subplot(n, n, i)
        x=[1,2]


        ### mf1 female
        y=[mf1_m.loc[genes[k],'GCKO2_d0'],mf1_m.loc[genes[k],'GCKO2_d13']]
        yerr=[mf1_sd.loc[genes[k],'GCKO2_d0'],mf1_sd.loc[genes[k],'GCKO2_d13']]
        plt.errorbar(x,y, yerr=yerr, fmt='b--')

        ### mf2 female
        y=[mf2_m.loc[genes[k],'GCKO2_d0'],mf2_m.loc[genes[k],'GCKO2_d13']]
        yerr=[mf2_sd.loc[genes[k],'GCKO2_d0'],mf2_sd.loc[genes[k],'GCKO2_d13']]
        plt.errorbar(x,y, yerr=yerr, fmt='r--')

        ### mf1 male
        y=[mf1_m.loc[genes[k],'145480_d0'],mf1_m.loc[genes[k],'145480_d13']]
        yerr=[mf1_sd.loc[genes[k],'145480_d0'],mf1_sd.loc[genes[k],'145480_d13']]
        plt.errorbar(x,y, yerr=yerr, fmt='b-')

        ### mf2 male
        y=[mf2_m.loc[genes[k],'145480_d0'],mf2_m.loc[genes[k],'145480_d13']]
        yerr=[mf2_sd.loc[genes[k],'145480_d0'],mf2_sd.loc[genes[k],'145480_d13']]
        plt.errorbar(x,y, yerr=yerr, fmt='r-')
        plt.ylabel('log2 relative fitness')
        plt.title(labels)
        plt.legend(('mf1', 'mf2'))
        plt.xticks([1, 1.5, 2],['day0', '', 'day13'],fontsize=15)

        plt.ylim(-18, 1)

        k=k+1
    pdf.savefig(fig)
    pdf.close()





def getNewIdfromPrevID(pgenes):
    prev_to_new=pickle.load(open('/Users/vikash/Documents/Projects/DBInfo/prevTonew_PBANKA.pickle','rb'))
    db_df=pd.read_csv('/Users/vikash/Documents/Projects/DBInfo/PBANKA_id_conversion.txt', sep='\t')
    db_df=db_df.fillna('NA')
    not_found=[]
    geneConv={}

    for item in pgenes:
        if item in prev_to_new.keys():
            tmp= db_df[db_df['Gene ID']==prev_to_new[item]].copy()
            if tmp.empty:
                not_found.append(item)
                geneConv[item]=item
            else:

                if tmp['Gene Name or Symbol'].to_list()[0]=='NA':
                    geneConv[item]=tmp['Product Description'].to_list()[0]+'|'+prev_to_new[item]
                else:
                    geneConv[item]=tmp['Product Description'].to_list()[0]+'|'+tmp['Gene Name or Symbol'].to_list()[0]+'|'+prev_to_new[item]
    print('not found genes', not_found)
    return geneConv


def makeBarcodeCSV():

    ### we are going to make a csv file for barcoe sample

    df=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/plasmogem_all_berghei_targeting_vectors_and_barcodes_jan20.csv")

    input=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Fertility_screen/DataFiles/input_vector.txt', sep='\t')
    tmp=df.loc[:,['gene_id', 'barcode'] ].copy()
    tmp=tmp.drop_duplicates()
    ### find unique barcodes

    uniq_barcodes=tmp['barcode'].unique().tolist()
    uniq_genes=tmp['gene_id'].unique().tolist()

    store_multiple=[]
    multiple_barcode=[]
    barcode_gene={}
    for gene in uniq_genes:
        xx=tmp[tmp['gene_id']==gene]
        if xx.shape[0]>1:
            multiple_barcode.append(gene)
            store_multiple.append(xx['barcode'].to_list())
        elif xx.shape[0]==1:
            barcode_gene[gene]=xx['barcode'].to_list()[0]
        else:
            print('not found')


    ## get barcode based on PlasmoDB ID



    see_genes=[]
    for item in input['PbGEM-ID'].to_list():
        if df[df['PlasmoGEM_ID']==item].empty:
            see_genes.append(item)

    ######### write files

    out=open("barcode_file.txt", 'w')
    out.write("gene\tbarcode\n")
    for k,v in barcode_gene.items():
        out.write("%s\t%s\n"%(k,v))
    for i in range(len(multiple_barcode)):
        out.write("%s\t%s\n"%(multiple_barcode[i],'|'.join(store_multiple[i])))

    ######### write files

    out=open("barcode_gene_file.csv", 'w')
    out.write("barcode,gene\n")
    for k in uniq_barcodes:
        ## get genes
        tmp=df[df['barcode']==k]
        plasmoGEM=set(tmp['PlasmoGEM_ID'].to_list())
        pbanka=set(tmp['gene_id'].to_list())

        gene='|'.join(plasmoGEM|pbanka)

        out.write("%s,%s\n"%(k,gene))


def getGroupMap(list1,list2):
    d={}
    for i,item in enumerate(list2):
        if item not in d.keys():
            d[item]=[]
        d[item].append(list1[i])
    return d


def plotHistogramonError(gut_values,gut_values_err,pcr_time_df):

    #### gut_time_df:  is mean values at day 13


    ### gut_time_df_err1: is the sd at day 13

    ####
    input_df=pcr_time_df[0] # this is the input dataframe at time day0
    input_values=getValues(input_df)

    pdf=plot_folder+'Dissection_error.svg'

    plotfitScatterDayNoLine(gut_values,gut_values_err,'log2(relative abundances)','Relative error',pdf,'Dissection_error')


    ### plot histogram
    pdf = plot_folder + 'Distribution_input(D0)_new.svg'
    plothistNew(input_values, '', 'Frequency', pdf, 50, [-17, -2])
    import pdb;pdb.set_trace()


def plothistNew(values,xlab,ylab,pdf,bins,xlim):
    ## get relative abundance
    data = np.array(values)
    # remove nans
    data = data[~np.isnan(data)]
    fig= plt.figure(figsize=(10,8))
# ax = plt.gca()
    kwargs = dict(histtype='step', alpha=1, normed=True, bins=bins,facecolor='blue', linewidth=2)

    plt.hist(data, **kwargs)
    # plt.grid(True)
    plt.xlim(xlim)
    plt.xlabel(xlab,fontsize=20)
    plt.ylabel(ylab,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig(pdf,format="svg")


def plotfitScatterDayNoLine(x,y,xlab,ylab,pdf,title):
    fig= plt.figure(figsize=(10,8))
    ax = plt.gca()
    xx=[]
    yy=[]
    pp=[]
    xpxp=[]
    y3=[]
    for i in range(len(x)):
        xx.append(np.array(x[i]))
        yy.append(np.array(y[i]))
        z=np.polyfit(x[i], y[i], 1)
        p=np.poly1d(z)
        mini=np.min(x[i])
        maxi=np.max(x[i])
        pp.append(p)
        xp = np.linspace(mini, maxi, 100)
#         print(xp,pp)
        xpxp.append(xp)

        # # exponnetial fit
        # popt, pcov = curve_fit(exponenial_func, xx[i], yy[i], p0=(1, 1e-3, 1))
        # print(popt)
        # y3.append(exponenial_func(xp, *popt))
    # lines=ax.plot(xx[0], yy[0],'gray.', xx[1], yy[1], 'balck.')
    # p1=ax.scatter(xx[0], yy[0], color='#999999', marker='.')
    p1=ax.scatter(xx[0], yy[0], color='#000000', marker='.')
    # p2=ax.scatter(xx[1], yy[1], color='#000000', marker='.')


    #lines=ax.plot(xx[0], yy[0], 'b.', xpxp[0], pp[0](xpxp[0]), 'b-',xx[1], yy[1], 'g.', xpxp[1], pp[1](xpxp[1]), 'g-',linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=20)
    plt.ylabel(ylab,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    # plt.rc('xtick',labelsize=)
    # plt.rc('ytick',labelsize=15)
    plt.legend([p1],['d13'],fontsize=16);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);
    plt.xlim([-17,-2])
    plt.savefig(pdf,format="svg")

    plt.show()
if __name__ == '__main__':

    #makeBarcodeCSV()
    stepByStep_barSeqAnalysis()
