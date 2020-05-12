import os
import sys

code=os.getcwd()
upLevel=code.replace('codes','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/Figures')

# print(sys.path)

from special_fun_fertility  import *
# input files which will be needed for analysis of cliare screen data

data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]

# ID conversion: we used plasmoDB to convert ID for P. Berghai
prev_to_new=pickle.load(open(data_folder+'/prevTonew_PBANKA.pickle','rb'))
db_df=pd.read_csv(data_folder+'/PBANKA_id_conversion.txt', sep='\t')
db_df=db_df.fillna('NA')
## end of databse information



def pre_process_pilot(manifests_df,count_df):

    ''' We will take average of the two reads '''
    # manifests_df.set_index('NGI Sample ID',inplace=True)
    req_df=count_df.copy()
    manifests_df=manifests_df.fillna('NA')

    manifests_df.set_index('Sample description',inplace=True)
    final_df=pd.DataFrame(index=req_df.index,columns=['Gene',  'Barcodes'])

    for ngi in manifests_df.index:

       ### get two reads
       # sam= manifests_df.loc[ngi,'Sample description'] # sample description
       sam=ngi
       two_reads=req_df.columns[req_df.columns.str.contains(sam)]
       if len(two_reads)==2:
           # final_df[ngi]=np.nan()
           final_df.loc[:,sam]=req_df.loc[:,two_reads].mean(axis=1)
           # final_df_two_read[ngi+'.read1']=np.nan()
           # final_df_two_read[ngi+'.read2']=np.nan()
           # final_df_two_read.loc[:,sam+'.read1']=req_df.loc[:,two_reads[0]]
           # final_df_two_read.loc[:,sam+'.read2']=req_df.loc[:,two_reads[1]]
       elif len(two_reads)==1:
           # final_df[ngi]=np.nan()
           final_df.loc[:,sam]=req_df.loc[:,two_reads].copy()
           print('Please check there is only one reads for the sample: %s'%sam)
       else:
           print('Number of reads are mismatched for the sample: %s'%sam)

    return final_df,manifests_df

def stepwiseAnalysis():
    ''' We are going do analyis of pilot study for Claire data'''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/data_pilot/manifest_pilot_small.txt",sep='\t')

    d28573=pd.read_csv(data_folder+"/data_pilot/counts_28573.csv")
    d28551=pd.read_csv(data_folder+"/data_pilot/counts_28551.csv")
    ### combine two data frames
    dfn = pd.merge(d28573, d28551, on='gene', how='outer')  # combined first and second files
    df_modi = dfn  #
    df_modi = df_modi.drop(df_modi.index[0])
    df_modi = df_modi.drop(['barcode_x', 'barcode_y'], axis=1)  # drop out some unnecessary columns
    df_modi.set_index(['gene'], inplace=True)
    count_df=df_modi.copy()
    input_df=pd.read_csv(data_folder+'/data_pilot/PbSTM168_cross_phenotypes_final.csv', sep=';')

    comm_genes=set(input_df['Gene_ID'])&set(count_df.index)

    filter_count_df=count_df.loc[comm_genes,:].copy()

    ### filter count df for pilot



    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_df,manfest_df=pre_process_pilot(manifests_df,filter_count_df)

    ## remove genes
    remove_genes=['PBANKA_051200']

    filtered_count_df = final_df.drop(remove_genes)
    filtered_count_df.to_csv(out_folder+"/filterd_count_matrix_pilot1.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    ## if you we do not want to plot then plot_info=None
    #plot_info={'pdf':out_folder+"/relative_abundance_of_pool1.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    #plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis
    
    error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)







    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)




if __name__ == '__main__':
    stepwiseAnalysis()
