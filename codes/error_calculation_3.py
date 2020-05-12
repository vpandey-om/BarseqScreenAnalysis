import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;
fig_save='/Users/vikash/Documents/Projects/Claire/Plots/'



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


def plothist(values,xlab,ylab,pdf,bins,xlim):
    ## get relative abundance
    data = np.array(values)
    # remove nans
    data = data[~np.isnan(data)]
    fig= plt.figure(figsize=(10,8))
# ax = plt.gca()
    plt.hist(data, bins)
    plt.xlim(xlim)
    plt.xlabel(xlab,fontsize=15)
    plt.ylabel(ylab,fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig(pdf)

def plothistNP(np_data, xlab, ylab, pdf, bins):
    ## get relative abundance
    data =np_data
    # remove nans
    data = data[~np.isnan(data)]
    fig = plt.figure(figsize=(10, 8))
    # ax = plt.gca()
    plt.hist(data, bins)
    plt.xlabel(xlab, fontsize=15)
    plt.ylabel(ylab, fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.savefig(pdf)

    # plt.show()

def getValues(df):
    df_values=[]
    for t_values in df.values:
        for val in t_values:
            df_values.append(val)
    return df_values



def readMeta():
    meta=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/meta2.csv')
    meta.set_index('Unnamed: 0', inplace=True)
    meta.index.names = ['Label']
    return meta

def cal_mean_and_sd_groupby_columns(df_modi,mapping):
    tmp_mean_T = pd.DataFrame(index=df_modi.index)
    tmp_std_T = pd.DataFrame(index=df_modi.index)
    for k,item in mapping.items():

        tmp_mean_T[k]=df_modi[item].mean(axis=1).copy()
        tmp_std_T[k]=df_modi[item].std(axis=1).copy()
    return tmp_mean_T,tmp_std_T

def cal_mean_and_sd_groupby_columns_old(df_modi):
    tmp_mean_T = pd.DataFrame(index=df_modi.index)
    tmp_std_T = pd.DataFrame(index=df_modi.index)
    for k, item in df_modi.groupby(df_modi.columns, axis=1).groups.items():
        tmp_mean_T[k]=df_modi[item].mean(axis=1).copy()
        tmp_std_T[k]=df_modi[item].std(axis=1).copy()
    return tmp_mean_T,tmp_std_T

def getGroupMap(list1,list2):
    d={}
    for i,item in enumerate(list2):
        if item not in d.keys():
            d[item]=[]
        d[item].append(list1[i])
    return d


def df_for_errors():
    d28573=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/counts_28573.csv")
    d28551=pd.read_csv("/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/counts_28551.csv")
    dfn = pd.merge(d28573, d28551, on='gene', how='outer')  # combined first and second files
    # import numpy as np
    df_modi = dfn  # This df is going to be modified
    # remove one row pandas for no match
    # df_modi.drop([0,1])
    df_modi = df_modi.drop(df_modi.index[0])
    df_modi = df_modi.drop(['barcode_x', 'barcode_y'], axis=1)  # drop out some unnecessary columns
    df_modi.set_index(['gene'], inplace=True)
    column_df=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/column_modi_df.txt','\t')
    df_modi.columns=column_df['Column2']
    df_ori=df_modi.copy()
    # calculate mean and satndared deviation for PCR replicates
    saveColumns = df_modi.columns

    pairedReadColumns = df_modi.columns.str.replace("\\.[12]$", '', regex=True)  # this is for pair sequencing ends
    df_modi.columns = pairedReadColumns

    # remove PbSTM168_
    PbSTM168 = df_modi.columns.str.replace("PbSTM168_", '', regex=True)
    df_modi.columns = PbSTM168
    map_paired_reads=getGroupMap(saveColumns, df_modi.columns)
    # paired_df = df_modi.copy() # this dataframe will we use to find diffrence between PCR

    # to find mf1.1 vas mf1.2 and  mf2.1 vs mf2.2
    pcr_columns = df_modi.columns.str.replace("_r[12]", '', regex=True)
    df_modi.columns = pcr_columns
    # pcr_df=df_modi.copy()
    map_pcr=getGroupMap(saveColumns, df_modi.columns)
    # find midgut columns
    mossifeed_df=df_modi.copy()
    midgut_flag = df_modi.columns.str.extract("^[0-9]+_(.+mf[12].[12])$")
    midgut_columns = saveColumns[np.where(~midgut_flag.isnull())[0]]
    # dissection error
    midgut_df = df_ori.loc[:, midgut_columns]

    midgut_df.columns = midgut_df.columns.str.replace("_r[12].[12]$", '', regex=True)
    midgut_df.columns = midgut_df.columns.str.replace(".[12]$", '', regex=True)

    # remove PbSTM168_
    PbSTM168 = midgut_df.columns.str.replace("PbSTM168_[0-9]+_", '', regex=True)
    midgut_df.columns = PbSTM168

    map_midgut = getGroupMap(midgut_columns, midgut_df.columns)
    # find mosquito feed diffrence mf1 and mf2

    mossifeed_flag = mossifeed_df.columns.str.extract("^[0-9]+_(.+mf[12])+")
    mossifeed_columns = saveColumns[np.where(~ mossifeed_flag.isnull())[0]]
    mossifeed_df = df_ori.loc[:, mossifeed_columns]

    #
    colNames1=mossifeed_df.columns.str.replace("\\.[12]$", '',regex=True)
    mossifeed_df.columns=colNames1
    colNames2=mossifeed_df.columns.str.replace("_r[12]", '',regex=True)
    mossifeed_df.columns=colNames2
    colNames4=mossifeed_df.columns.str.replace("\\.[12]$", '', regex=True)
    mossifeed_df.columns=colNames4
    colNames5=mossifeed_df.columns.str.replace("\\_mf[12]$", '', regex=True)
    mossifeed_df.columns=colNames5
    colNames3=mossifeed_df.columns.str.replace("PbSTM168_[0-9]+_", '', regex=True)
    mossifeed_df.columns=colNames3

    # colNames4=mossifeed_df.columns.str.replace("^[0-9]+_(.+mf[12]).+$", "\\1",regex=True)
    # mossifeed_df.columns=colNames4

    # mossifeed_df.columns = mossifeed_df.columns.str.replace("_r[12].[12]$", '', regex=True)
    #
    # # remove PbSTM168_
    # PbSTM168 = mossifeed_df.columns.str.replace("PbSTM168_", '', regex=True)
    # mossifeed_df.columns = PbSTM168
    # mossifeed_df.columns=mossifeed_df.columns.str.replace('[^A-Za-z0-9_]+[12]', '', regex=True)
    map_mossifeed = getGroupMap(mossifeed_columns, mossifeed_df.columns)

    return df_ori,map_pcr,map_midgut,map_mossifeed


def calRel(df):
    n=df.shape[1]
    df=df.div(df.sum(axis=1), axis=0)
    # get log value
    # # cut=1e-10 # this is the cutoff
    # df[df<cut]=cut
    df=np.log2(df)
    #df=np.log10(df)
    return df


def calRel_withoutLog(df):

    n=df.shape[1]
    df=df.div(df.sum(axis=1), axis=0)
    # get log value
    # cut=1e-10
    # df[df<cut]=cut;
    return df

def plotfitScatter(x,y,xlab,ylab,pdf,title):
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
    lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], pp[0](xpxp[0]), 'r-',xx[1], yy[1], 'b.', xpxp[1], pp[1](xpxp[1]), 'b-',xx[2], yy[2], 'g.', xpxp[2], pp[2](xpxp[2]), 'g-',linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=15)
    plt.ylabel(ylab,fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(lines,['d0','d0','d7','d7','d14','d14'],fontsize=10);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);
#     plt.xlim([-14,1])
    plt.savefig(pdf)

    plt.show()


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
    lines=ax.plot(xx[0], yy[0], 'b.', xpxp[0], pp[0](xpxp[0]), 'b-',xx[1], yy[1], 'g.', xpxp[1], pp[1](xpxp[1]), 'g-',linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=15)
    plt.ylabel(ylab,fontsize=15)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(lines,['d7','d7','d14','d14'],fontsize=10);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);

    plt.savefig(pdf)

    plt.show()



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
    p1=ax.scatter(xx[0], yy[0], color='#999999', marker='.')
    p2=ax.scatter(xx[1], yy[1], color='#000000', marker='.')


    #lines=ax.plot(xx[0], yy[0], 'b.', xpxp[0], pp[0](xpxp[0]), 'b-',xx[1], yy[1], 'g.', xpxp[1], pp[1](xpxp[1]), 'g-',linewidth=3)
    # lines=ax.plot(xx[0], yy[0], 'r.', xpxp[0], y3[0], 'r-',xx[1], yy[1], 'b.', xpxp[1], y3[1], 'b-',xx[2], yy[2], 'g.', xpxp[2], y3[2], 'g-')
#     lines=ax.plot( xpxp[0], pp[0](xpxp[0]), 'r-', xpxp[1], pp[1](xpxp[1]), 'b-', xpxp[2], pp[2](xpxp[2]), 'g-')
    plt.xlabel(xlab,fontsize=20)
    plt.ylabel(ylab,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    # plt.rc('xtick',labelsize=)
    # plt.rc('ytick',labelsize=15)
    plt.legend([p1,p2],['d7','d14'],fontsize=16);
    plt.title(title)
#     plt.legend(lines,['d0','d7','d14'],fontsize=15);
    plt.xlim([-17,-2])
    plt.savefig(pdf,format="svg")

    plt.show()


def calErrorEachTime(df_m, df_sd,mapping,time,conds):
    #with this function we want to compute error for each time point
    # df_m is the mean df_sd is SD.
    # test whether there is duplicates in genes
    if len(df_m.index) !=len(df_m.index.unique()):
        print ("there is duplicate genes: We should fixed which gene do we want to take in analysis")
        exit()

    comm=readMutPhenotypes(df_m) # common mutants

    df_m=df_m.loc[comm,:] # this is just to take only mutants
    df_sd=df_sd.loc[comm,:] # this is just to take only mutants
    mapOrigIdx=findMatch(df_m.columns,mapping)
    # import pdb;pdb.set_trace()

    time_df=[]
    time_df_err=[]
    # technical replicates

    for t in time:
        tech_df = pd.DataFrame(index=df_m.index)
        tech_df_err = pd.DataFrame(index=df_sd.index)
        for cond in conds:
            tech_df = pd.concat([tech_df,df_m[mapOrigIdx[(cond,t)]].copy()],axis=1)
            tech_df_err = pd.concat([tech_df_err,df_sd[mapOrigIdx[(cond,t)]].copy()],axis=1) #  error
        cut=1e-10
        if tech_df.isnull().values.any():
            print ("relative abundance contains nan value do you want to ignore? we put 1")
            tech_df=tech_df.fillna(1)
        if tech_df_err.isnull().values.any():
            print ("relative abundance error contains nan value do you want to ignore? we put 0")
            tech_df_err=tech_df_err.fillna(0)
        # import pdb;pdb.set_trace()
        tech_df[tech_df < cut]=cut
        tech_df=np.log2(tech_df)
        tech_df_err[tech_df_err < cut]=cut
        tech_df_err=np.log2(tech_df_err)
        time_df.append(tech_df)
        time_df_err.append(tech_df_err)



    # collaspe in one vector calculate errors and values
    l_values=[]
    for item in time_df:
# technical replicates ratio
        values=getValues(item)
        l_values.append(values)

    l_values_err=[]
    for item in time_df_err:
# technical replicates error
        values=getValues(item)
        l_values_err.append(values)

    return time_df,time_df_err,l_values,l_values_err


def getColumnsFormDF(df,tag):
    t_cols=[]
    for t in tag:
        flag = df.columns.str.contains(t)
        cols= df.columns[np.where(flag)[0]]
        t_cols.append(cols)
    return t_cols

    # a column



def selectColumns(df_m, df_sd,tag):
    #with this function we want to compute error for each time point
    # df_m is the mean df_sd is SD.
    # test whether there is duplicates in genes
    if len(df_m.index) !=len(df_m.index.unique()):
        print ("there is duplicate genes: We should fixed which gene do we want to take in analysis")
        exit()

    comm=readMutPhenotypes(df_m) # common mutants

    df_m=df_m.loc[comm,:] # this is just to take only mutants
    df_sd=df_sd.loc[comm,:] # this is just to take only mutants
    df_rel=df_sd/df_m # this is just to take only mutants
    time_cols=getColumnsFormDF(df_m,tag)

    # import pdb;pdb.set_trace()
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

def findMatch(main_list,mapping):

    mapOrigIdx={}
    for k,item in mapping.items():
        index_list=[]
        for subs in item:
            res = list(filter(lambda x: subs in x, main_list))
            for idx in res:
                index_list.append(idx)
        mapOrigIdx[k]=index_list
    return mapOrigIdx

def plotAllDistributions(ref_ori,df_ori):
    pdf=fig_save+'allRatio_distribution.pdf'
    vals=getValues(ref_ori)
    plothist(vals,'ratios','frequency',pdf,500)

    pdf=fig_save+'allcounts_distribution.pdf'
    vals=getValues(df_ori)
    plothist(vals,'counts','frequency',pdf,500)

def readMutPhenotypes(df):
    anot=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/PbSTM168_cross_phenotypes_final.csv',';')
    s1=set(list(df.index))
    #anot["newGene"] = anot["Gene_Nme"].map(str) +'_' + anot["E_phenotype"]
    anot["newGene"] = anot["Gene_Nme"].map(str)
    # anot.head()
    d1=pd.Series(anot['newGene'].values,index=anot['Gene_ID']).to_dict()
    l1=[]
    l2=[]
    d2={}
    for k,item in d1.items():
        l1.append(k.strip())
        l2.append(item.strip())
        d2[k.strip()]=item.strip()
    comm=s1&set(l1)
    return comm



def calRelAbundance():
    # calulate mean and standared deviation for three types of error
    df_ori,map_pcr,map_midgut,map_mossifeed=df_for_errors()


    # oliver suggested put count 1 if it is less than one
    # df_ori=df_ori.fillna(0)
    # df_ori[df_ori<0.5]=0.5

    df_ori=df_ori+1
    # import pdb;pdb.set_trace()
    ref_ori=calRel_withoutLog(df_ori.T)
    ref_ori=ref_ori.T


    # plot distributions
    #plotAllDistributions(ref_ori,df_ori)

    # hey we get same distributions


    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr) # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed)  # this is the mean and difference between mosiquito feeds

    meta=readMeta()
    mossiTimemap={} # this is the dictionary for mosquito feed time points
    for k,item in meta.groupby(['parasite_background','timepoint']).indices.items():
        mossiTimemap[k]=meta.index[item]
    # calculate relative abundance for pcr
    conds=['Cl15cy1','145480_GOMO','GCKO2_GOMO']

    #calculate error for PCR duplicates
    time=['mf0','mf7','mf14']
    pcr_time_df,pcr_time_df_err,pcr_values,pcr_values_err=calErrorEachTime(pcr_m, pcr_sd ,mossiTimemap,time,conds) # find only for mutants
    pdf=fig_save+'PCR_error_relative abundance.pdf'

    plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    # between the midGuts

    time=['mf7','mf14']
    gut_time_df,gut_time_df_err,gut_values,gut_values_err=calErrorEachTime( midgut_m, midgut_sd ,mossiTimemap,time,conds) # find only for mutants
    pdf=fig_save+'midgut_error_relative abundance.pdf'
    plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    import pdb;pdb.set_trace()
    time=['mf0','mf7','mf14']
    mossifeed_time_df,mossifeed_time_df_err,mossifeed_values,mossifeed_values_err=calErrorEachTime( mossifeed_m, mossifeed_sd ,mossiTimemap,time,conds) # find only for mutants
    pdf=fig_save+'mossifeed_error_relative abundance.pdf'
    plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Mosquito feed error')




def calRelAbundanceRedo():
    # calulate mean and standared deviation for three types of error
    df_ori,map_pcr,map_midgut,map_mossifeed=df_for_errors()

    comm = readMutPhenotypes(df_ori)  # common mutants

    df_ori=df_ori.loc[comm,:]
    # import pdb;pdb.set_trace()
    # oliver suggested put count 1 if it is less than one
    # df_ori=df_ori.fillna(0)
    # df_ori[df_ori<1]=1
    df_ori=df_ori+1
    # df_ori=df_ori
    # import pdb;pdb.set_trace()
    ref_ori=calRel_withoutLog(df_ori.T)
    ref_ori=ref_ori.T

    # plot distributions
    #plotAllDistributions(ref_ori,df_ori)

    # hey we get same distributions
    import pdb;pdb.set_trace()

    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr) # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed)  # this is the mean and difference between mosiquito feeds

    #calculate error for PCR duplicates
    tag=['d0_mf','d7_mf','d14_mf']

    pcr_time_df,pcr_time_df_err,pcr_time_df_err1,pcr_time_df_nl,pcr_values,pcr_values_err,pcr_values_err1,pcr_values_nl=selectColumns(pcr_m, pcr_sd,tag)
    pdf=fig_save+'PCR_error_relative abundance.pdf'

    #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','SD of relative abundances',pdf,'PCR_error')
    pdf=fig_save+'PCR_error_relative abundance vs rel variance.pdf'
    # import pdb;pdb.set_trace()
    #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    plotfitScatter(pcr_values,pcr_values_err1,'log2(mean of relative abundances)','(SD of relative abundances)/mean',pdf,'PCR_error')

    ####
    tag=['d7_mf','d14_mf']
    gut_time_df,gut_time_df_err,gut_time_df_err1,gut_time_df_nl,gut_values,gut_values_err,gut_values_err1,gut_values_nl=selectColumns( midgut_m, midgut_sd,tag) # find only for mutants

    pdf=fig_save+'midgut_error_relative abundance.pdf'
    #plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','SD of relative abundances',pdf,'Dissection_error')
    pdf=fig_save+'disection_error.pdf'
    pdf=fig_save+'disection_error.svg'
    import pdb;pdb.set_trace()
    #plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    # plotfitScatterDay(gut_values,gut_values_err1,'log2(mean relative abundances)','SD of relative abundances/Mean',pdf,'Dissection_error')
    plotfitScatterDayNoLine(gut_values,gut_values_err1,'log2(relative abundances)','Relative error',pdf,'Dissection_error')

    ####
    tag=['d0','d7','d14']
    mossifeed_time_df,mossifeed_time_df_err,mossifeed_time_df_err1,mossifeed_time_df_nl,mossifeed_values,mossifeed_values_err,mossifeed_values_err1,mossifeed_values1=selectColumns( mossifeed_m, mossifeed_sd ,tag) # find only for mutants

    pdf=fig_save+'mossifeed_error_relative abundance.pdf'
    #plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Mosquito feed error')
    plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','SD of relative abundances',pdf,'Mosquito feed error')

    pdf=fig_save+'mossifeed_error_relative abundance vs rel varaince.pdf'
    #plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Mosquito feed error')
    plotfitScatter(mossifeed_values,mossifeed_values_err1,'log2(mean relative abundances)','SD of relative abundances/Mean',pdf,'Mosquito feed error')


def testDist():
        # we want to test distribution of count and relative frequency
    df_ori, map_pcr, map_midgut, map_mossifeed = df_for_errors()

    comm = readMutPhenotypes(df_ori)  # common mutants
    df_ori = df_ori.loc[comm, :]




    ########## Plot count distribution #####
    #df_p = df_ori
    # df_p.columns = df_p.columns.str.replace('PbSTM168_', '')
    # pairedReadColumns = df_p.columns.str.replace("\\.[12]$", '', regex=True)  # this is for pair sequencing ends
    # df_p.columns = pairedReadColumns
    # map_read = getGroupMap(df_ori.columns, pairedReadColumns)
    # dfp_m, dfp_s = cal_mean_and_sd_groupby_columns(df_ori, map_read)
    # dfp_m = dfp_m + 0.01
    # dfp_m = np.log2(dfp_m)
    # fig=dfp_m.hist(figsize=(50, 50))
    # [x.title.set_size(12) for x in fig.ravel()]
    # plt.savefig(fig_save+'distibution_columns.pdf')
    # L=getValues(df_ori )
    # plt.hist(L, 20)

    # plt.show()

    df_ori=df_ori+1

        # take distribution of input and plot on disection error

    ref_ori = calRel_withoutLog(df_ori.T)
    ref_ori = ref_ori.T

        # plot distributions
        # plotAllDistributions(ref_ori,df_ori)

        # hey we get same distributions


    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr)  # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed)  # this is the mean and difference between mosiquito feeds

    # calculate error for PCR duplicates
    tag = ['d0_mf', 'd7_mf', 'd14_mf']

    pcr_time_df, pcr_time_df_err, pcr_time_df_err1, pcr_time_df_nl,pcr_values, pcr_values_err, pcr_values_err1,pcr_values_nl = selectColumns(pcr_m, pcr_sd, tag)


    ####
    tag = ['d7_mf', 'd14_mf']
    gut_time_df, gut_time_df_err, gut_time_df_err1, gut_time_df_nl,gut_values, gut_values_err, gut_values_err1,gut_values_nl = selectColumns(midgut_m, midgut_sd, tag)  # find only for mutants


    joinplot(gut_values, gut_values_err1, pcr_values,pcr_values_nl)


def joinplot(gut_values, gut_values_err1,pcr_values,pcr_values_nl):
    # get min and max axis
    aa = np.array(gut_values)
    maxi=aa.max()
    mini=aa.min()
    aa=np.array(pcr_values[0])
    max1=aa.max()
    min1=aa.min()
    if min1<mini:
        mini=min1
    if max1>maxi:
        maxi=max1
    aa = np.array(pcr_values_nl[0])
    n = len(pcr_values_nl[0])
    mean_val=aa.mean()
    sd_val = aa.std()

    # mean_val=np.log2(mean_val)
    # sd_val = np.log2(sd_val)


    #pdf = fig_save + 'Distribution_Day0.pdf'
    pdf = fig_save + 'Distribution_Day0_new.svg'
    plothistNew(pcr_values[0], '', 'Frequency', pdf, 50, [-17, -2])

    size=2*n
    pdf = fig_save + 'Distribution_Day0_size2.pdf'

    L = np.random.normal(mean_val, sd_val, size)
    ln=np.log2(L)

    # plothistNP(ln, 'relative abundance (Day 0)', 'Distribution (size %d)' % size, pdf, 20)

    p2 = L / np.sum(L)
    p2 = np.log2(p2)
    pdf = fig_save + 'Distribution_Day0_size2_p2.pdf'
    # plothistNP(p2, 'relative abundance (Day 0)', 'Distribution (size %d)' % size, pdf, 20)

    pdf = fig_save + 'Distribution_Day0_size2_p1.pdf'

    # L = np.array(pcr_values_nl[0])
    L = np.random.normal(mean_val, sd_val, n)
    p1 = L / np.sum(L)
    p1 = np.log2(p1)

    # plothistNP(p1, 'relative abundance (Day 0)', 'Distribution (size %d)' % n, pdf, 20)

    size2=3*n
    L = np.random.normal(mean_val, sd_val, size2)
    p3 = L / np.sum(L)
    p3 = np.log2(p3)


    size3 = 4 * n
    L = np.random.normal(mean_val, sd_val, size3)
    p4 = L / np.sum(L)
    p4 = np.log2(p4)

    ######
    bins=30
    plt.hist(p1, bins, alpha=0.5, label='size=%d'%n)
    plt.hist(p2, bins, alpha=0.5, label='2*size=%d'%size)
    plt.hist(p3, bins, alpha=0.5, label='3*size=%d' %size2)
    plt.hist(p4, bins, alpha=0.5, label='4*size=%d' % size3)
    plt.legend(loc='upper right')
    plt.show()

def weightedMeanAnalysis():
    # here we compute weighted mean and variance
      # we want to test distribution of count and relative frequency
    df_ori, map_pcr, map_midgut, map_mossifeed = df_for_errors()
    comm = readMutPhenotypes(df_ori)  # common mutants
    df_ori = df_ori.loc[comm, :]
    df_ori = df_ori + 1

    # take distribution of input and plot on disection error

    ref_ori = calRel_withoutLog(df_ori.T)
    ref_ori = ref_ori.T
    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr)  # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    cmb_fitness,cmb_fitness_without_log,tl=calmidgutRatio(pcr_m, pcr_sd)


    ## we should plot scatater plot
    plotrelFitness(cmb_fitness_without_log)


def calmidgutRatio(df_m,df_sd):
    # based on pcr mean and sd we will calculate ratios
    # df_sd = df_sd / df_m  # this is just to take only mutants
    time= ['d0', 'd7', 'd14']
    mf=['mf1','mf2']
    backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']
    ratioDict={}
    for b in backgrounds:
        for t in time:
            for f in mf:
                string=b+'_'+t+'_'+f
                time_cols = getColumnsFormDF(df_m, [string])

                ratioDict[(b,t,f)]=time_cols[0]

    time = ['d7', 'd14']
    rel_fit=pd.DataFrame(index=df_m.index)
    rel_fit1 = pd.DataFrame(index=df_m.index)
    rel_err = pd.DataFrame(index=df_m.index)
    rel_err1 = pd.DataFrame(index=df_m.index) #without log scale
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
                rel_fit1 = pd.concat([rel_fit1, tech_df.iloc[:, 0:].div(input.iloc[:, 0], axis=0)], axis=1)  # error

                # calculate error of A/B
                term1=tech_df.iloc[:, 0:].div(input.iloc[:, 0], axis=0).pow(2,axis=1)
                term2= tech_df_err.div(tech_df,axis=1).pow(2, axis=1)
                term3=input_err.div(input,axis=1).pow(2, axis=1)
                term4=term2.iloc[:,0:].add(term3.iloc[:, 0],axis=0)
                final_term=term1.iloc[:,0:].multiply(term4.iloc[:,0],axis=0)

                rel_err1 = pd.concat([rel_err1, final_term], axis=1)  # error
                # calculate relative frequency
                tech_df = np.log2(tech_df)
                input = np.log2(input)

                rel_fit = pd.concat([rel_fit, tech_df.iloc[:, 0:].sub(input.iloc[:, 0], axis=0)], axis=1)  # error
                # SD on log2 scale
                tech_df_err=np.log2(tech_df_err)
                input_err=np.log2(input_err)

                tmp=tech_df_err.pow(2, axis=1).iloc[:,0:].add(input_err.pow(2,axis=1).iloc[:,0],axis=0)
                # rel_err=pd.concat([rel_err,tmp.apply(np.sqrt)],axis=1) SD
                rel_err=pd.concat([rel_err,tmp],axis=1) # varaiance


    ## distributions of relative fitness


    dfp_m = np.log2(rel_fit1)
    fig=dfp_m.hist(figsize=(50, 50))
    [x.title.set_size(12) for x in fig.ravel()]
    # plt.savefig(fig_save+'distibution_of_ratios.pdf')
    #
    # import pdb;pdb.set_trace()
    cmb_fitness={}
    cmb_fitness_without_log={}
    #
    # for b in backgrounds:
    #     for t in time:
    #         string=b+'_'+t
    #         tcols = getColumnsFormDF(rel_fit, [string])
    #         tm=rel_fit[tcols[0]].copy()
    #         terr=rel_err[tcols[0]].copy()
    #         cmb_fitness[(b,t)]=gaussianMeanAndVariance(tm,terr)
    #
    #         ## without log
    #
    #         tcols = getColumnsFormDF(rel_fit1, [string])
    #         tm=rel_fit1[tcols[0]].copy()
    #         terr=rel_err1[tcols[0]].copy()
    #         cmb_fitness_without_log[(b,t)]=gaussianMeanAndVariance(tm,terr)

    return cmb_fitness,cmb_fitness_without_log,[(rel_fit, rel_err),(rel_fit1, rel_err1)]


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





def getshapeColor(index,ess_dict, sex_dict):

    import matplotlib.pyplot as plt
    markers=[]
    for i in index:
        markers.append("o")
    # get markers ess_dict
    for k,item in ess_dict.items():
        if k=='essential':
            for i in item:
                markers[i]="D"
        elif k=='non-essential':
            for i in item:
                markers[i]="s"
        else:
            continue
    m_colors=[]
    # get color
    for i in index:
        m_colors.append("gray")
    for k,item in sex_dict.items():
        if k=='male':
            for i in item:
                m_colors[i]="blue"
        elif k=='female':
            for i in item:
                m_colors[i]="orange"
        else:
            continue
    return m_colors,markers

def plotrelFitness(cmb_fitness):
    anot=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/PbSTM168_cross_phenotypes_final.csv',';')
    anot["newGene"] = anot["Gene_Nme"].map(str)
    d1=pd.Series(anot['newGene'].values,index=anot['Gene_ID']).to_dict()

    anot.set_index('Gene_ID',inplace=True)
    comm=set(cmb_fitness[('Cl15cy1STM', 'd14')][0].keys()) & set(anot.index)
    t_anot=anot.loc[comm,:]
    #### plot

    ess_dict=t_anot.groupby('Essential_phenotype').indices
    sex_dict=t_anot.groupby('Published_cross_phenotype').indices
    m_colors,markers=getshapeColor(comm,ess_dict, sex_dict)

    labels=[]
    for item in comm:
        if item in d1.keys():
            labels.append(d1[item])
        else:
            labels.append(item)

    #backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']
    ## take in comm order

    normal=[]
    male=[]
    female=[]

    # cmb_fitness[('Cl15cy1STM', 'd14')][0]=np.log10(cmb_fitness[('Cl15cy1STM', 'd14')][0])
    # cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0]=np.log10(cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0])
    # cmb_fitness[('145480GOMOSTM', 'd14')][0]=np.log10(cmb_fitness[('145480GOMOSTM', 'd14')][0])

    for idx in comm:
        normal.append(cmb_fitness[('Cl15cy1STM', 'd14')][0][idx])
        male.append(cmb_fitness[('145480GOMOSTM', 'd14')][0][idx])
        female.append(cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0][idx])


    import pdb;pdb.set_trace()

    pdf=fig_save+'145480_GOMO_GCKO2_GOMO_test2.pdf'
    plot2dScatter(male, female,labels,'145480_GOMO','GCKO2_GOMO',pdf,markers,m_colors)

    pdf=fig_save+'Cl15cy1_GCKO2_GOMO_test2.pdf'
    plot2dScatter(normal, female,labels,'Cl15cy1','GCKO2_GOMO',pdf,markers,m_colors)

    pdf=fig_save+'Cl15cy1_145480_GOMO_test2.pdf'
    plot2dScatter(male, normal,labels,'145480_GOMO','Cl15cy1',pdf,markers,m_colors)


def plot2dScatter(fx, fy,labels,x_lab,y_lab,pdf_path,markers,m_colors):
    fig= plt.figure(figsize=(15,15))
    ax = plt.gca()


#plot errorbars
    for i in np.arange(0, len(fx)):
        ax.plot(fx[i], fy[i], linestyle="None", marker=markers[i],ms=15,color=m_colors[i])
        # ax.plot([fx[i]+xerror[i], fx[i]-xerror[i]], [fy[i], fy[i]], marker="_",color='black', linewidth=2)
        # ax.plot([fx[i], fx[i]], [fy[i]+yerror[i], fy[i]-yerror[i]], marker="_",color='black', linewidth=2)

    cnt=0
    for i,j in zip(fx, fy):
        corr = -0.010 # adds a little correction to put annotation in marker's centrum
#         ax.annotate(labels[cnt],  xy=(i + corr, j + corr),fontsize=12)
        ax.text(i + corr, j + corr,labels[cnt],
        verticalalignment='bottom', horizontalalignment='right',
        color='black', fontsize=10)
        cnt=cnt+1
    #configure axes
#     ax.legend()
#     plt.xlim(-5, 1)
#     plt.ylim(-5, 1)
    plt.xlabel(x_lab,fontsize=25)
    plt.ylabel(y_lab,fontsize=25)
#     plt.legend(['o','s','v'], ['NA', 'NE', 'E'],loc='upper left',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.savefig(pdf_path)
    plt.show()



def calDissectionRatio(df_m,df_sd):
    # based on pcr mean and sd we will calculate ratios
    # df_sd = df_sd / df_m  # this is just to take only mutants
    time= ['d0', 'd7', 'd14']
    mf=['mf1','mf2']
    backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']
    ratioDict={}
    for b in backgrounds:
        for t in time:
            for f in mf:
                string=b+'_'+t+'_'+f
                time_cols = getColumnsFormDF(df_m, [string])

                ratioDict[(b,t,f)]=time_cols[0]

    time = ['d7', 'd14']
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
                # rel_fit1 = pd.concat([rel_fit1, tech_df.iloc[:, 0:].div(input.iloc[:, 0], axis=0)], axis=1)  # error
                #
                # # calculate error of A/B
                # term1=tech_df.iloc[:, 0:].div(input.iloc[:, 0], axis=0).pow(2,axis=1)
                # term2= tech_df_err.div(tech_df,axis=1).pow(2, axis=1)
                # term3=input_err.div(input,axis=1).pow(2, axis=1)
                # term4=term2.iloc[:,0:].add(term3.iloc[:, 0],axis=0)
                # final_term=term1.iloc[:,0:].multiply(term4.iloc[:,0],axis=0)
                #
                # rel_err1 = pd.concat([rel_err1, final_term], axis=1)  # error
                # calculate relative frequency
                # tech_df = np.log2(tech_df)
                # input = np.log2(input)

                rel_fit = pd.concat([rel_fit, tech_df.iloc[:, 0:].sub(input.iloc[:, 0], axis=0)], axis=1)  # error

                # SD on log2 scale
                # tech_df_err=np.log2(tech_df_err)
                # input_err=np.log2(input_err)

                tmp=tech_df_err.pow(2, axis=1).iloc[:,0:].add(input_err.pow(2,axis=1).iloc[:,0],axis=0)
                # rel_err=pd.concat([rel_err,tmp.apply(np.sqrt)],axis=1) SD
                rel_err=pd.concat([rel_err,tmp],axis=1) # varaiance


    ## distributions of relative fitness


    # dfp_m = np.log2(rel_fit1)
    # fig=dfp_m.hist(figsize=(50, 50))
    # [x.title.set_size(12) for x in fig.ravel()]
    # plt.savefig(fig_save+'distibution_of_ratios.pdf')
    #
    # import pdb;pdb.set_trace()

    #
    # for b in backgrounds:
    #     for t in time:
    #         string=b+'_'+t
    #         tcols = getColumnsFormDF(rel_fit, [string])
    #         tm=rel_fit[tcols[0]].copy()
    #         terr=rel_err[tcols[0]].copy()
    #         cmb_fitness[(b,t)]=gaussianMeanAndVariance(tm,terr)
    #
    #         ## without log
    #
    #         tcols = getColumnsFormDF(rel_fit1, [string])
    #         tm=rel_fit1[tcols[0]].copy()
    #         terr=rel_err1[tcols[0]].copy()
    #         cmb_fitness_without_log[(b,t)]=gaussianMeanAndVariance(tm,terr)

    return rel_fit, rel_err



def step_wise_test():
    df_ori, map_pcr, map_midgut, map_mossifeed = df_for_errors()
    comm = readMutPhenotypes(df_ori)  # common mutants
    df_ori = df_ori.loc[comm, :]
    df_ori = df_ori + 1

    # take distribution of input and plot on disection error

    # ref_ori = calRel_withoutLog(df_ori.T)
    ref_ori = calRel(df_ori.T) # with log2 scale

    ref_ori = ref_ori.T
    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr)  # this is the mean and difference between pcr duplicates

    # control genes
    # cmb_fitness,cmb_fitness_without_log,tl=calmidgutRatio(pcr_m, pcr_sd)
    rel_fitness,rel_var= calDissectionRatio(pcr_m, pcr_sd)

    # withlogtuple=tl[0]
    # nologtuple=tl[1]
    # normalized by reference poposed by oliver
    control_genes=['PBANKA_102460','PBANKA_100210','PBANKA_100220','PBANKA_050120']


    # calculate  inverse weighted mean for control genes
    # rel_fitness=withlogtuple[0]
    # rel_var=withlogtuple[1]



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
    backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']
    time = ['d7', 'd14']
    for b in backgrounds:
        for t in time:
            string=b+'_'+t
            tcols = getColumnsFormDF(normalized_fit, [string])
            tm=normalized_fit[tcols[0]].copy()
            terr=normalized_var[tcols[0]].copy()

            cmb_fitness[(b,t)]=gaussianMeanAndVariance(tm,terr)

    #plotrelFitness(cmb_fitness)
    ##########  This function  is used for preparing plots for fertility website
    res_df=prepareForPGfertility(cmb_fitness)
    import pdb;pdb.set_trace

    ############################
    viz_df=calPvalZscore(cmb_fitness)

    # plotInSns(viz_df)
    plotInSnsNew(viz_df)

def readAnotation():
    anot=pd.read_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/PbSTM168_cross_phenotypes_final.csv',';')
    anot["newGene"] = anot["Gene_Nme"].map(str)
    d1 = pd.Series(anot['newGene'].values, index=anot['Gene_ID']).to_dict()
    anot.set_index('Gene_ID', inplace=True)

    return anot,d1


def calPvalZscore(cmb_fitness):
    import scipy.stats as st
    import scipy.special as sp
    anot, geneMap=readAnotation()
    # create df for ploting
    for k in cmb_fitness.keys():
        key1=k
        break


    comm=set(cmb_fitness[key1][0].index)& set(anot.index)
    geneName=[geneMap[item] for item in comm]
    viz_df=anot.loc[comm,:].copy()
    # viz_df['geneName'] = geneName

    pval=0.001
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


def prepareForPGfertility(cmb_fitness):
    # we are going to prepare for fertility
    import scipy.stats as st  # Statistical library
    import scipy.special as sp
    anot, geneMap=readAnotation()
    # create df for ploting
    for k in cmb_fitness.keys():
        key1=k
        break


    comm=set(cmb_fitness[key1][0].index)& set(anot.index)
    geneName=[geneMap[item] for item in comm]
    viz_df=anot.loc[comm,:].copy()
    # viz_df['geneName'] = geneName

    pval=0.001
    pval1=0.01
    upcut=np.log2(1)
    lowcut=np.log2(0.03)

    pheno_pval={}

    for key,item in cmb_fitness.items():

    # we will calculate p-vlaue and z score
        m=item[0] #  this is the relative fitness
        s=item[1] # varaiance
        z=(upcut-m)/s  #  z score
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

    viz_df['145480GOMOSTM_d14_pheno']=viz_df['145480GOMOSTM_d14_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
    viz_df['GCKO#2GOMOSTM_d14_pheno']=viz_df['GCKO#2GOMOSTM_d14_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
    viz_df.to_csv('/Users/vikash/Documents/Projects/Claire/Pilot_fertility_screen_data/Fertility_data_vis.csv',sep='\t')
    import pdb;pdb.set_trace()  # we are checking the pvalues for normal genes .
def plotInSns(viz_df):

    sns.set(style="whitegrid")

    print('hi')

    # ['Gene_Nme', 'Essential_phenotype', 'Predicted_cross_phenotype',
    #  'Published_cross_phenotype', 'newGene', 'Cl15cy1STM_d7_pheno',
    #  'Cl15cy1STM_d7_rel', 'Cl15cy1STM_d14_pheno', 'Cl15cy1STM_d14_rel',
    #  'GCKO#2GOMOSTM_d7_pheno', 'GCKO#2GOMOSTM_d7_rel',
    #  'GCKO#2GOMOSTM_d14_pheno', 'GCKO#2GOMOSTM_d14_rel',
    #  '145480GOMOSTM_d7_pheno', '145480GOMOSTM_d7_rel',
    #  '145480GOMOSTM_d14_pheno', '145480GOMOSTM_d14_rel'],

    viz_df.columns=viz_df.columns.str.replace('STM','')
    viz_df.columns = viz_df.columns.str.replace('rel', 'RGR')
    # change df_values

    viz_df['145480GOMO_d14_pheno']=viz_df['145480GOMO_d14_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
    viz_df['GCKO#2GOMO_d14_pheno']=viz_df['GCKO#2GOMO_d14_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
    viz_df['Published_cross_phenotype']=viz_df['Published_cross_phenotype'].replace({'N': 'NA'})
    viz_df['Cl15cy1_d14_pheno']=viz_df['Cl15cy1_d14_pheno'].replace({'E': 'I', 'NE': 'F', 'NA': 'RF'})

    # for d7

    viz_df['145480GOMO_d7_pheno'] = viz_df['145480GOMO_d7_pheno'].replace({'E': 'IM', 'NE': 'FM', 'NA': 'RM'})
    viz_df['GCKO#2GOMO_d7_pheno'] = viz_df['GCKO#2GOMO_d7_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})

    viz_df['Cl15cy1_d7_pheno'] = viz_df['Cl15cy1_d7_pheno'].replace({'E': 'I', 'NE': 'F', 'NA': 'RF'})


    ###################### BAR PLOTS BEGINS

    # ## bar plots
    # bar_df=viz_df.sort_values(by=["145480GOMO_d14_RGR"])
    # plt.figure(figsize=(13, 20))
    # # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", hue="145480GOMO_d14_pheno", data=viz_df,ci=None,palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))
    # # plt.errorbar(x=viz_df["145480GOMO_d14_RGR"], y=ax.get_yticks(),xerr=viz_df["145480GOMO_d14_sd"], fmt='none', c='gray')
    # cmap={'FM':"#66c2a5", 'RM':"#8da0cb", 'IM':"#fc8d62"}
    # clrs = [cmap[x]  for x in bar_df["145480GOMO_d14_pheno"]]
    #
    # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", data=bar_df, ci=None,palette = clrs)
    # plt.errorbar(x=bar_df["145480GOMO_d14_RGR"], y=ax.get_yticks(), xerr=bar_df["145480GOMO_d14_sd"], fmt='none',c='gray')
    #
    #
    # # plt.legend(["#66c2a5","#8da0cb","#fc8d62"], ['FM', 'RM', 'IM'])
    # fig = ax.get_figure()
    # fig.savefig(fig_save + "bar_male_plot_RGR_.pdf")
    # # fig, ax = plt.subplots()
    # # plot = viz_df["145480GOMO_d14_RGR"].plot(kind='bar', yerr=viz_df["145480GOMO_d14_sd"], colormap='OrRd_r', edgecolor='black', grid=False, figsize=(8, 8),
    # #                     ax=ax, position=0.45, error_kw=dict(ecolor='black', elinewidth=0.5), width=0.8)
    # # plt.show()
    #
    #
    #
    #
    # ## bar plots
    # bar_df = viz_df.sort_values(by=["GCKO#2GOMO_d14_RGR"])
    # plt.figure(figsize=(13, 20))
    # # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", hue="145480GOMO_d14_pheno", data=viz_df,ci=None,palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))
    # # plt.errorbar(x=viz_df["145480GOMO_d14_RGR"], y=ax.get_yticks(),xerr=viz_df["145480GOMO_d14_sd"], fmt='none', c='gray')
    # cmap = {'FF': "#66c2a5", 'RF': "#8da0cb", 'IF': "#fc8d62"}
    # clrs = [cmap[x] for x in bar_df["GCKO#2GOMO_d14_pheno"]]
    #
    # ax = sns.barplot(x="GCKO#2GOMO_d14_RGR", y="newGene", data=bar_df, ci=None, palette=clrs)
    # plt.errorbar(x=bar_df["GCKO#2GOMO_d14_RGR"], y=ax.get_yticks(), xerr=bar_df["GCKO#2GOMO_d14_sd"], fmt='none',c='gray')
    #
    # # plt.legend(["#66c2a5","#8da0cb","#fc8d62"], ['FM', 'RM', 'IM'])
    # fig = ax.get_figure()
    # fig.savefig(fig_save + "bar_female_plot_RGR_.pdf")


     ###################### BAR PLOTS END

    remove_genes=['PBANKA_143240','PBANKA_051200']
    import pdb;pdb.set_trace()
    # remove indexes
    tmp_df=viz_df.copy()
    tmp_df=tmp_df.drop(remove_genes)
    plt.figure(figsize=(8, 8))
    # viz_df.groupby(["145480GOMOSTM_d14_pheno","GCKO#2GOMOSTM_d14_pheno"]).indices

    ax= sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR",hue = "145480GOMO_d14_pheno",markers = dict(FF="o", RF="s", IF="X"), style = "GCKO#2GOMO_d14_pheno", data = viz_df,legend='brief',\
    palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))

    mapLabel={"145480GOMO_d14_pheno":'Male mutant',"GCKO#2GOMO_d14_pheno":'Female mutant'}
    orders=['145480GOMO_d14_pheno', 'FM', 'RM', 'IM', 'GCKO#2GOMO_d14_pheno', 'FF', 'RF', 'IF']
    file=fig_save + "male_female_pval_based_RGR_.pdf"
    set_handles_labels(ax,mapLabel,orders,viz_df,file,'145480GOMO_d14_RGR','GCKO#2GOMO_d14_RGR','Male mutant','Female mutant')
    # ax.legend(loc='lower left',fancybox=True, framealpha=0.3)
    # label_point(viz_df['145480GOMO_d14_RGR'], viz_df['GCKO#2GOMO_d14_RGR'], viz_df['newGene'], plt.gca())
    # fig = ax.get_figure()
    # fig.savefig(fig_save+"pval_based_RGR.pdf")


    # g = sns.FacetGrid(viz_df, col="145480GOMOSTM_d14_pheno", hue="GCKO#2GOMOSTM_d14_pheno")
    # g = (g.map(plt.scatter, "145480GOMOSTM_d14_rel", "GCKO#2GOMOSTM_d14_rel", edgecolor="w").add_legend())
    # plt.show()
    plt.figure(figsize=(8, 8))
    # ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype",
    #                      style="Essential_phenotype", data=viz_df, legend='brief')
    ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype", data=viz_df,legend='brief', \
                         palette=dict(female="#66c2a5", NA="#8da0cb", male="#fc8d62"))
    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df['145480GOMO_d14_RGR'], viz_df['GCKO#2GOMO_d14_RGR'], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "published_RGR.pdf")

    # plt.legend(loc='upper left')
    # plt.show()

    ## visualize based on known phenotypes

    # male vs wild type

    plt.figure(figsize=(8, 8))
    ax = sns.scatterplot(x="Cl15cy1_d14_RGR", y="145480GOMO_d14_RGR", hue="145480GOMO_d14_pheno",markers=dict(F="o", RF="s", I="X"), style="Cl15cy1_d14_pheno", data=viz_df,legend='brief', \
                         palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))

    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df["Cl15cy1_d14_RGR"], viz_df["145480GOMO_d14_RGR"], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "Cl15cy1_vs_male_pval_based_RGR_.pdf")

    # female vs wild type
    plt.figure(figsize=(8, 8))
    ax = sns.scatterplot(x="GCKO#2GOMO_d14_RGR", y = "Cl15cy1_d14_RGR", hue = "GCKO#2GOMO_d14_pheno", markers = dict(F="o", RF="s",I="X"), style = "Cl15cy1_d14_pheno", data = viz_df, legend = 'brief', \
    palette = dict(FF="#66c2a5", RF="#8da0cb", IF="#fc8d62"))

    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df["GCKO#2GOMO_d14_RGR"], viz_df["Cl15cy1_d14_RGR"], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "Cl15cy1_vs_female_pval_based_RGR_.pdf")


def plotInSnsNew(viz_df):
    # import pdb;pdb.set_trace()

    # sns.set_style("white")
    sns.set(style="white",font_scale=1.5)
    rm_dict=removeGenes()
    print('hi')
    # ['GCKO#2', '145480', 'Cl15cy1']
    # ['Gene_Nme', 'Essential_phenotype', 'Predicted_cross_phenotype',
    #  'Published_cross_phenotype', 'newGene', 'Cl15cy1STM_d7_pheno',
    #  'Cl15cy1STM_d7_rel', 'Cl15cy1STM_d14_pheno', 'Cl15cy1STM_d14_rel',
    #  'GCKO#2GOMOSTM_d7_pheno', 'GCKO#2GOMOSTM_d7_rel',
    #  'GCKO#2GOMOSTM_d14_pheno', 'GCKO#2GOMOSTM_d14_rel',
    #  '145480GOMOSTM_d7_pheno', '145480GOMOSTM_d7_rel',
    #  '145480GOMOSTM_d14_pheno', '145480GOMOSTM_d14_rel'],

    viz_df.columns=viz_df.columns.str.replace('STM','')
    viz_df.columns = viz_df.columns.str.replace('rel', 'RGR')
    # change df_values

    viz_df['145480GOMO_d14_pheno']=viz_df['145480GOMO_d14_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
    viz_df['GCKO#2GOMO_d14_pheno']=viz_df['GCKO#2GOMO_d14_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
    viz_df['Published_cross_phenotype']=viz_df['Published_cross_phenotype'].replace({'N': 'NA'})
    viz_df['Cl15cy1_d14_pheno']=viz_df['Cl15cy1_d14_pheno'].replace({'E': 'I', 'NE': 'F', 'NA': 'RF'})

    # for d7

    viz_df['145480GOMO_d7_pheno'] = viz_df['145480GOMO_d7_pheno'].replace({'E': 'IM', 'NE': 'FM', 'NA': 'RM'})
    viz_df['GCKO#2GOMO_d7_pheno'] = viz_df['GCKO#2GOMO_d7_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})

    viz_df['Cl15cy1_d7_pheno'] = viz_df['Cl15cy1_d7_pheno'].replace({'E': 'I', 'NE': 'F', 'NA': 'RF'})


    # ##################### BAR PLOTS BEGINS
    #
    # ## bar plots
    # bar_df=viz_df.sort_values(by=["145480GOMO_d14_RGR"])
    #
    # # remove genes
    # bar_df=bar_df.drop(rm_dict['145480'])
    #
    # plt.figure(figsize=(16, 20))
    # # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", hue="145480GOMO_d14_pheno", data=viz_df,ci=None,palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))
    # # plt.errorbar(x=viz_df["145480GOMO_d14_RGR"], y=ax.get_yticks(),xerr=viz_df["145480GOMO_d14_sd"], fmt='none', c='gray')
    # cmap={'FM':"#66c2a5", 'RM':"#8da0cb", 'IM':"#fc8d62"}
    # cmap={'FM':"#1b9e77", 'RM':"#7570b3", 'IM':"#d95f02"}
    # clrs = [cmap[x]  for x in bar_df["145480GOMO_d14_pheno"]]
    #
    # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", data=bar_df, ci=None,palette = clrs)
    # plt.errorbar(x=bar_df["145480GOMO_d14_RGR"], y=ax.get_yticks(), xerr=bar_df["145480GOMO_d14_sd"], fmt='none',c='gray')
    #
    #
    # # plt.legend(["#66c2a5","#8da0cb","#fc8d62"], ['FM', 'RM', 'IM'])
    # fig = ax.get_figure()
    # fig.savefig(fig_save + "bar_145480_plot_RGR.pdf")
    # # fig, ax = plt.subplots()
    # # plot = viz_df["145480GOMO_d14_RGR"].plot(kind='bar', yerr=viz_df["145480GOMO_d14_sd"], colormap='OrRd_r', edgecolor='black', grid=False, figsize=(8, 8),
    # #                     ax=ax, position=0.45, error_kw=dict(ecolor='black', elinewidth=0.5), width=0.8)
    # # plt.show()
    #
    #
    #
    #
    # ## bar plots
    # bar_df = viz_df.sort_values(by=["GCKO#2GOMO_d14_RGR"])
    #  # remove genes
    # bar_df=bar_df.drop(rm_dict['GCKO#2'])
    #
    # plt.figure(figsize=(16, 20))
    # # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", hue="145480GOMO_d14_pheno", data=viz_df,ci=None,palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))
    # # plt.errorbar(x=viz_df["145480GOMO_d14_RGR"], y=ax.get_yticks(),xerr=viz_df["145480GOMO_d14_sd"], fmt='none', c='gray')
    # cmap={'FF':"#1b9e77", 'RF':"#7570b3", 'IF':"#d95f02"}
    # clrs = [cmap[x] for x in bar_df["GCKO#2GOMO_d14_pheno"]]
    #
    # ax = sns.barplot(x="GCKO#2GOMO_d14_RGR", y="newGene", data=bar_df, ci=None, palette=clrs)
    # plt.errorbar(x=bar_df["GCKO#2GOMO_d14_RGR"], y=ax.get_yticks(), xerr=bar_df["GCKO#2GOMO_d14_sd"], fmt='none',c='gray')
    #
    # # plt.legend(["#66c2a5","#8da0cb","#fc8d62"], ['FM', 'RM', 'IM'])
    # fig = ax.get_figure()
    # fig.savefig(fig_save + "bar_GCKO_plot_RGR_.pdf")
    #
    # plt.figure(figsize=(16, 20))
    # ## for normal background
    # bar_df = viz_df.sort_values(by=["Cl15cy1_d14_RGR"])
    #  # remove genes
    # bar_df=bar_df.drop(rm_dict['Cl15cy1'])
    # # ax = sns.barplot(x="145480GOMO_d14_RGR", y="newGene", hue="145480GOMO_d14_pheno", data=viz_df,ci=None,palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))
    # # plt.errorbar(x=viz_df["145480GOMO_d14_RGR"], y=ax.get_yticks(),xerr=viz_df["145480GOMO_d14_sd"], fmt='none', c='gray')
    # cmap={'F':"#1b9e77", 'RF':"#7570b3", 'I':"#d95f02"}
    # clrs = [cmap[x] for x in bar_df["Cl15cy1_d14_pheno"]]
    #
    # ax = sns.barplot(x="Cl15cy1_d14_RGR", y="newGene", data=bar_df, ci=None, palette=clrs)
    # plt.errorbar(x=bar_df["Cl15cy1_d14_RGR"], y=ax.get_yticks(), xerr=bar_df["Cl15cy1_d14_sd"], fmt='none',c='gray')
    #
    # # plt.legend(["#66c2a5","#8da0cb","#fc8d62"], ['FM', 'RM', 'IM'])
    # fig = ax.get_figure()
    # fig.savefig(fig_save + "bar_Cl15cy1_plot_RGR_.pdf")
    #
    #
    #  ###################### BAR PLOTS END
    import pdb;pdb.set_trace()
    remove_genes=[]

    for k,v in rm_dict.items():
        if k in ['GCKO#2', '145480']:
            for it in v:
                if it not in remove_genes:
                    remove_genes.append(it)
    # remove indexes
    tmp_df=viz_df.copy()
    tmp_df=tmp_df.drop(remove_genes)
    NAlist=tmp_df.index[tmp_df['Published_cross_phenotype'] == 'NA'].tolist()
    tmp_df['newGene'][NAlist]=''

    ## plot published one
    plt.figure(figsize=(8, 8))
    # ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype",
    #                      style="Essential_phenotype", data=viz_df, legend='brief')
    # ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype", data=viz_df,legend='brief', \
    #                      palette=dict(female="#66c2a5", NA="#8da0cb", male="#fc8d62"))

    # ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", style="Published_cross_phenotype", data=tmp_df,legend='brief', \
    #                      markers = dict(male="o", female='^', NA="X"), s=60)
    # plt.legend(loc='lower right')

    ##################

    # Set style of scatterplot
    sns.set_context("notebook", font_scale=1.1)
    sns.set_style("ticks")

    # Create scatterplot of dataframe
    sns.lmplot(x="145480GOMO_d14_RGR", # Horizontal axis
           y="GCKO#2GOMO_d14_RGR", # Vertical axis
           data=tmp_df, # Data source
           fit_reg=False, # Don't fix a regression line
           hue="Published_cross_phenotype", # Set color
           scatter_kws={"marker": "D", # Set marker style
                        "s": 100}) # S marker size

# Set title
    plt.title('Histogram of IQ')

    # Set x-axis label
    plt.xlabel('Time')

    # Set y-axis label
    plt.ylabel('Deaths')








    #################
    import pdb;pdb.set_trace()

    # file=fig_save + "male_female_published_RGR_no_legend.pdf"
    # fig = ax.get_figure()
    # fig.savefig(file)
    # # ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    # label_point(viz_df['145480GOMO_d14_RGR'], viz_df['GCKO#2GOMO_d14_RGR'], viz_df['newGene'], plt.gca())
    mapLabel={"Published_cross_phenotype":'Phenotype','NA':'No data','female':'Female','male':'Male'}
    orders=['Published_cross_phenotype', 'female', 'male', 'NA']
    file=fig_save + "male_female_published_RGR.pdf"

    set_handles_labels(ax,mapLabel,orders,tmp_df,file,'145480GOMO_d14_RGR','GCKO#2GOMO_d14_RGR','log2(female KO relative growth rate)','log2(Male KO relative growth rate)')





    ##

    plt.figure(figsize=(8, 8))
    # viz_df.groupby(["145480GOMOSTM_d14_pheno","GCKO#2GOMOSTM_d14_pheno"]).indices

    ax= sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR",hue = "145480GOMO_d14_pheno",markers = dict(FF="o", RF="s", IF="X"), style = "GCKO#2GOMO_d14_pheno", data = viz_df,legend='brief',\
    palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))

    mapLabel={"145480GOMO_d14_pheno":'Male mutant',"GCKO#2GOMO_d14_pheno":'Female mutant'}
    orders=['145480GOMO_d14_pheno', 'FM', 'RM', 'IM', 'GCKO#2GOMO_d14_pheno', 'FF', 'RF', 'IF']
    file=fig_save + "male_female_pval_based_RGR_.pdf"
    set_handles_labels(ax,mapLabel,orders,viz_df,file,'145480GOMO_d14_RGR','GCKO#2GOMO_d14_RGR','Male mutant','Female mutant')
    # ax.legend(loc='lower left',fancybox=True, framealpha=0.3)
    # label_point(viz_df['145480GOMO_d14_RGR'], viz_df['GCKO#2GOMO_d14_RGR'], viz_df['newGene'], plt.gca())
    # fig = ax.get_figure()
    # fig.savefig(fig_save+"pval_based_RGR.pdf")


    # g = sns.FacetGrid(viz_df, col="145480GOMOSTM_d14_pheno", hue="GCKO#2GOMOSTM_d14_pheno")
    # g = (g.map(plt.scatter, "145480GOMOSTM_d14_rel", "GCKO#2GOMOSTM_d14_rel", edgecolor="w").add_legend())
    # plt.show()
    plt.figure(figsize=(8, 8))
    # ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype",
    #                      style="Essential_phenotype", data=viz_df, legend='brief')
    ax = sns.scatterplot(x="145480GOMO_d14_RGR", y="GCKO#2GOMO_d14_RGR", hue="Published_cross_phenotype", data=viz_df,legend='brief', \
                         palette=dict(female="#66c2a5", NA="#8da0cb", male="#fc8d62"))
    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df['145480GOMO_d14_RGR'], viz_df['GCKO#2GOMO_d14_RGR'], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "published_RGR.pdf")

    # plt.legend(loc='upper left')
    # plt.show()

    ## visualize based on known phenotypes

    # male vs wild type

    plt.figure(figsize=(8, 8))
    ax = sns.scatterplot(x="Cl15cy1_d14_RGR", y="145480GOMO_d14_RGR", hue="145480GOMO_d14_pheno",markers=dict(F="o", RF="s", I="X"), style="Cl15cy1_d14_pheno", data=viz_df,legend='brief', \
                         palette=dict(FM="#66c2a5", RM="#8da0cb", IM="#fc8d62"))

    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df["Cl15cy1_d14_RGR"], viz_df["145480GOMO_d14_RGR"], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "Cl15cy1_vs_male_pval_based_RGR_.pdf")

    # female vs wild type
    plt.figure(figsize=(8, 8))
    ax = sns.scatterplot(x="GCKO#2GOMO_d14_RGR", y = "Cl15cy1_d14_RGR", hue = "GCKO#2GOMO_d14_pheno", markers = dict(F="o", RF="s",I="X"), style = "Cl15cy1_d14_pheno", data = viz_df, legend=Flase, \
    palette = dict(FF="#66c2a5", RF="#8da0cb", IF="#fc8d62"))

    ax.legend(loc='lower left', fancybox=True, framealpha=0.3)
    label_point(viz_df["GCKO#2GOMO_d14_RGR"], viz_df["Cl15cy1_d14_RGR"], viz_df['newGene'], plt.gca())
    fig = ax.get_figure()
    fig.savefig(fig_save + "Cl15cy1_vs_female_pval_based_RGR_.pdf")



def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.06, point['y'], str(point['val']),fontsize=6)



def plotBasedOnPvale(cmb_fitness):
    anot = pd.read_csv('/Users/vikash/Documents/Claire/Codes/PbSTM168_cross_phenotypes_final.csv', ';')
    anot["newGene"] = anot["Gene_Nme"].map(str)
    d1=pd.Series(anot['newGene'].values,index=anot['Gene_ID']).to_dict()

    anot.set_index('Gene_ID',inplace=True)
    comm=set(cmb_fitness[('Cl15cy1STM', 'd14')][0].keys()) & set(anot.index)
    t_anot=anot.loc[comm,:]
    #### plot

    ess_dict=t_anot.groupby('Essential_phenotype').indices
    sex_dict=t_anot.groupby('Published_cross_phenotype').indices
    m_colors,markers=getshapeColor(comm,ess_dict, sex_dict)

    labels=[]
    for item in comm:
        if item in d1.keys():
            labels.append(d1[item])
        else:
            labels.append(item)

    #backgrounds=['Cl15cy1STM','GCKO#2GOMOSTM','145480GOMOSTM']
    ## take in comm order

    normal=[]
    male=[]
    female=[]

    # cmb_fitness[('Cl15cy1STM', 'd14')][0]=np.log10(cmb_fitness[('Cl15cy1STM', 'd14')][0])
    # cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0]=np.log10(cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0])
    # cmb_fitness[('145480GOMOSTM', 'd14')][0]=np.log10(cmb_fitness[('145480GOMOSTM', 'd14')][0])

    for idx in comm:
        normal.append(cmb_fitness[('Cl15cy1STM', 'd14')][0][idx])
        male.append(cmb_fitness[('145480GOMOSTM', 'd14')][0][idx])
        female.append(cmb_fitness[('GCKO#2GOMOSTM', 'd14')][0][idx])


    import pdb;pdb.set_trace()

    pdf=fig_save+'145480_GOMO_GCKO2_GOMO_test2.pdf'
    plot2dScatter(male, female,labels,'145480_GOMO','GCKO2_GOMO',pdf,markers,m_colors)

    pdf=fig_save+'Cl15cy1_GCKO2_GOMO_test2.pdf'
    plot2dScatter(normal, female,labels,'Cl15cy1','GCKO2_GOMO',pdf,markers,m_colors)

    pdf=fig_save+'Cl15cy1_145480_GOMO_test2.pdf'
    plot2dScatter(male, normal,labels,'145480_GOMO','Cl15cy1',pdf,markers,m_colors)


def set_handles_labels(ax,mapLabel,orders,viz_df,file,df_xlab,df_ylab,xlab,ylab):
    handles,labels=ax.get_legend_handles_labels()
    new_handles=[]
    new_labels=[]
    for o in orders:
        idx=labels.index(o)
        new_handles.append(handles[idx])
        if o in mapLabel.keys():
            new_labels.append(mapLabel[o])
        else:
            new_labels.append(o)
    # ax.legend(loc='upper right', fancybox=True, framealpha=0.1)

    label_point(viz_df[df_xlab], viz_df[df_ylab], viz_df['newGene'], plt.gca())
    # plt.legend(bbox_to_anchor=(1.1, 1.05), loc='lower left', ncol=1)
    plt.legend(loc='lower left')
    ax.legend(labels=new_labels,handles=new_handles)

    import pdb;pdb.set_trace()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    # ax.set_xlim([-13, 2.5])
    # ax.set_ylim([-13, 2.5])
    fig = ax.get_figure()
    fig.savefig(file)

    print('plotting....')

def overlay():
    from PIL import Image

    background = Image.open("/Users/vikash/Documents/Projects/Claire/Plots/Distribution_Day0.svg")
    overlay = Image.open("/Users/vikash/Documents/Projects/Claire/Plots/disection_error.svg")

    background = background.convert("RGBA")
    overlay = overlay.convert("RGBA")

    new_img = Image.blend(background, overlay, 0.5)
    new_img.save("/Users/vikash/Documents/Projects/Claire/Plots/new.png","PNG")




def calErrorAndDist():
    # calulate mean and standared deviation for three types of error
    df_ori,map_pcr,map_midgut,map_mossifeed=df_for_errors()

    comm = readMutPhenotypes(df_ori)  # common mutants

    df_ori=df_ori.loc[comm,:]
    # import pdb;pdb.set_trace()
    # oliver suggested put count 1 if it is less than one
    # df_ori=df_ori.fillna(0)
    # df_ori[df_ori<1]=1
    df_ori=df_ori+1
    # df_ori=df_ori
    # import pdb;pdb.set_trace()
    ref_ori=calRel_withoutLog(df_ori.T)
    ref_ori=ref_ori.T

    pcr_m, pcr_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_pcr) # this is the mean and difference between pcr duplicates

    midgut_m, midgut_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_midgut)  # this is the mean and difference between midugut samples

    # mossifeed_m, mossifeed_sd = cal_mean_and_sd_groupby_columns(ref_ori,map_mossifeed)  # this is the mean and difference between mosiquito feeds

    # #calculate error for PCR duplicates
    tag=['d0_mf','d7_mf','d14_mf']

    pcr_time_df,pcr_time_df_err,pcr_time_df_err1,pcr_time_df_nl,pcr_values,pcr_values_err,pcr_values_err1,pcr_values_nl=selectColumns(pcr_m, pcr_sd,tag)
    # import pdb;pdb.set_trace()
    # pdf=fig_save+'PCR_error_relative abundance.pdf'

    # #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    # plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','SD of relative abundances',pdf,'PCR_error')
    # pdf=fig_save+'PCR_error_relative abundance vs rel variance.pdf'
    # import pdb;pdb.set_trace()
    #plotfitScatter(pcr_values,pcr_values_err,'log2(mean of relative abundances)','log2(SD of relative abundances)',pdf,'PCR_error')
    # plotfitScatter(pcr_values,pcr_values_err1,'log2(mean of relative abundances)','(SD of relative abundances)/mean',pdf,'PCR_error')

    ####
    tag=['d7_mf','d14_mf']
    gut_time_df,gut_time_df_err,gut_time_df_err1,gut_time_df_nl,gut_values,gut_values_err,gut_values_err1,gut_values_nl=selectColumns( midgut_m, midgut_sd,tag) # find only for mutants

    # pdf=fig_save+'midgut_error_relative abundance.pdf'
    # #plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    # plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','SD of relative abundances',pdf,'Dissection_error')
    # pdf=fig_save+'disection_error.pdf'
    # pdf=fig_save+'disection_error.svg'
    #
    # #plotfitScatterDay(gut_values,gut_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Dissection_error')
    # # plotfitScatterDay(gut_values,gut_values_err1,'log2(mean relative abundances)','SD of relative abundances/Mean',pdf,'Dissection_error')
    # plotfitScatterDayNoLine(gut_values,gut_values_err1,'log2(relative abundances)','Relative error',pdf,'Dissection_error')
    #
    # ####
    # tag=['d0','d7','d14']
    # mossifeed_time_df,mossifeed_time_df_err,mossifeed_time_df_err1,mossifeed_time_df_nl,mossifeed_values,mossifeed_values_err,mossifeed_values_err1,mossifeed_values1=selectColumns( mossifeed_m, mossifeed_sd ,tag) # find only for mutants
    #
    # pdf=fig_save+'mossifeed_error_relative abundance.pdf'
    # #plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Mosquito feed error')
    # plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','SD of relative abundances',pdf,'Mosquito feed error')
    #
    # pdf=fig_save+'mossifeed_error_relative abundance vs rel varaince.pdf'
    # #plotfitScatter(mossifeed_values,mossifeed_values_err,'log2(mean relative abundances)','log2(SD of relative abundances)',pdf,'Mosquito feed error')
    # plotfitScatter(mossifeed_values,mossifeed_values_err1,'log2(mean relative abundances)','SD of relative abundances/Mean',pdf,'Mosquito feed error')
    #
    return gut_time_df,gut_time_df_err1,pcr_time_df



def plotErrorAndDit():
    rm_dict=removeGenes() # remove genes
    gut_time_df,gut_time_df_err1,pcr_time_df=calErrorAndDist()
    index=[0,1] # 0 is day 7 1 is day 14
    days=['d7', 'd14']


    all_background_df=[]
    all_background_df_err=[]
    for ind in index:
        background_df=[]
        background_df_err=[]
        df=gut_time_df[ind].copy()
        df1=gut_time_df_err1[ind].copy()
        for k in rm_dict.keys():
            back_cols=getColumnsFormDF(gut_time_df[ind],[k])
            back_cols=back_cols[0]
            tmp=df.loc[:,back_cols].copy()
            tmp=tmp.drop(rm_dict[k])
            background_df.append(tmp)
            # for errors
            tmp=df1.loc[:,back_cols].copy()
            tmp=tmp.drop(rm_dict[k])
            background_df_err.append(tmp)
        all_background_df.append(background_df)
        all_background_df_err.append( background_df_err)

    # collect d7 and d14

    values=[]
    values_err=[]
    for ind in index:
        tmp=[]
        tmp_err=[]
        for b in range(3):
            vals=getValues(all_background_df[ind][b])

            for itm in vals:
                tmp.append(itm)
            vals=getValues(all_background_df_err[ind][b])
            for itm in vals:
                tmp_err.append(itm)

        values.append(tmp)
        values_err.append(tmp_err)

    # plot disection error

    pdf=fig_save+'Dissection_error.svg'

    plotfitScatterDayNoLine(values,values_err,'log2(relative abundances)','Relative error',pdf,'Dissection_error')

    # plot distribution
    # remove genes
    input_df=pcr_time_df[0] # this is the input dataframe at time day0
    input_values=[]
    for k in rm_dict.keys():
        back_cols=getColumnsFormDF(input_df,[k])
        back_cols=back_cols[0]
        tmp=input_df.loc[:,back_cols].copy()
        tmp=tmp.drop(rm_dict[k])
        vals=getValues(tmp)
        for it in vals:
            input_values.append(it)

    pdf = fig_save + 'Distribution_input(D0)_new.svg'
    plothistNew(input_values, '', 'Frequency', pdf, 50, [-17, -2])


def removeGenes():
    rm_df=pd.read_csv('/Users/vikash/Documents/Projects/Claire/low_input_genes_to_remove.txt',sep='\t')
    rm_dict={} # thiss is the dictionary for removing genes from background

    for col in rm_df.columns:
        rm_dict[col]= rm_df[col].dropna()
    return rm_dict




def regenFigures():
    # we are regenerating figures
    plotErrorAndDit()






if __name__ == "__main__":
    #calRelAbundanceRedo() #This script is to cal culate relative abundances
    #testDist()
    #weightedMeanAnalysis()
    step_wise_test()
    #regenFigures()
