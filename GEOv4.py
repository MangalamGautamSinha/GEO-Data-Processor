#pip install GEOparse pandas requests functools
import os
import GEOparse
import pandas as pd
import numpy as np
import requests
from functools import reduce

dfs=[] #list of all the data frames that are merged later (GSM and GPL)

def get_file(GSE_No):
    try:
        os.environ['GEOPARSE_USE_HTTP_FOR_FTP'] = 'yes'
        gse = GEOparse.get_GEO(geo=GSE_No, destdir="./Data")
        return gse
    except OSError as e:
        print(f"An error occurred: {e}")
        print('Please manually download '+GSE_No+'_family.soft.gz file and place it inside the Data folder (without extracting).')
        quit()

def trim_unigene(string):
    index = string.find('//')
    if index != -1:
        string = string[:index]
    return string if string!="---" else np.nan


def get_gene_symbol(input_list):
        input_string = ','.join(input_list) #converting the GB_Acc list to a single string
        headers = {'content-type': 'application/x-www-form-urlencoded'} # searching the GB_Acc string on mygene.info and storing the output gene symbol in gene_symbols_list
        params = 'q='+input_string+'&scopes=accession,unigene&species=human&fields=symbol&size=1'
        res = requests.post('http://mygene.info/v3/query', data=params, headers=headers)
        gene_symbol_list = []
        if res.status_code == 200:
            data=res.json()
            for entry in data:
                if 'symbol' in entry:
                    gene_symbol = entry['symbol']
                else:
                    gene_symbol = '--'
                gene_symbol_list.append(gene_symbol)
        else:
            print("Error:", res.status_code)
        return gene_symbol_list


def process_gpl(gse):
    for gpl in gse.gpls.values(): 
        print(gpl.table.columns[1:].tolist())
        add_col=input("Enter the name of the column to keep from the above list: ")
        gpl.table = gpl.table[['ID', add_col]] #removing all columns except ID and additional column from GPL table , then appending this dataframe to dfs list
        gpl.table['ID'] = gpl.table['ID'].astype(str)
        dfs.append(gpl.table)

        df=gpl.table.copy()        
        print('Adding the column and searching for Gene Symbol...')
        if add_col=="unigene":
            df['unigene']=df['unigene'].map(trim_unigene)            
            sub_df = df[~df['unigene'].isna()].copy() #creating a sub data frame which only contains ID and unigene, provided unigene is not empty. To search for less number of gene symbol
            unigene_list=sub_df['unigene'].to_list()
            sub_df.loc[:, 'Symbol'] = get_gene_symbol(unigene_list)
            dfs.append(sub_df)

            
        elif add_col=="GB_ACC":
            sub_df = df[~df['GB_ACC'].isna()].copy()
            gb_acc_list=sub_df['GB_ACC'].to_list()
            sub_df['Symbol'] = get_gene_symbol(gb_acc_list)
            dfs.append(sub_df.drop('GB_ACC', axis=1)) #dropping GB_ACC column to avoid duplicate column in final_df, as gpl.table already contain GB_ACC


            
def process_gsm(gse):
    k=0
    i=str
    for gsm_name, gsm in gse.gsms.items(): #appending all the GSM dataframes containing ID and GSM No. into dfs list
        if k==0:
            print(gsm.get_metadata_attribute("characteristics_ch1"))
            i=input('Enter the index (position from beginning) of subtype starting from 0, For Non TNBC enter N: ')
            k=1
        
        if i == 'N' or i == 'n':
            subtype=''
        else:
            subtype=gsm.get_metadata_attribute("characteristics_ch1")[int(i)]
            index = subtype.find(':')
            subtype = ' ('+subtype[index+2:]+')'
                              
        
        gsm.table.rename(columns={'ID_REF':'ID', 'VALUE':gsm_name+subtype}, inplace=True) #Rename column heading and add subtype in column name
        gsm.table['ID'] = gsm.table['ID'].astype(str)
        dfs.append(gsm.table)



def merge_dataframes(dfs):
    final_df = reduce(lambda  left,right: pd.merge(left,right,on=['ID'],how='outer'), dfs)
    last_col=final_df.columns[-1] #Dropping the whole row, if data at the last column is empty
    final_df.drop(final_df[final_df[last_col].isna()].index, inplace=True)
    return final_df


GSE_No=input("Enter GSE No.: ")
GSE_No=GSE_No if GSE_No[0:3]=='GSE' else 'GSE'+GSE_No
gse=get_file(GSE_No)

print("Processing GPLs")
process_gpl(gse)
print("Processing GSMs")
process_gsm(gse)
print("Merging Dataframes")
final_df=merge_dataframes(dfs)
final_df.to_csv(GSE_No+'.csv', index=False)
print("Completed")