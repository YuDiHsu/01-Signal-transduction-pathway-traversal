import pandas as pd
import os
import mygene
import pickle
from pprint import pprint


df_1 = pd.read_excel(os.path.join('.', 'imported_data', '四組都有數值.xlsx'),
                     index_col=None, header=1, sheet_name='4組都有數值')
df_2 = pd.read_excel(os.path.join('.', 'imported_data', 'Gene_expression_level_A.xlsx'),
                     index_col=None, header=1, sheet_name='FPKM&TPM')
df_1.loc[:, 'GN'] = df_1.loc[:, 'GN'].apply(lambda x: x.upper())
df_2.loc[:, 'gene_name'] = df_2.loc[:, 'gene_name'].apply(lambda x: x.upper())

df_1.loc[:, 'verify'] = 1

df_1 = df_1.loc[(df_1['nor Heavy / Light_pho-total'] <= -1) |
                (df_1['nor Heavy / Light_pho-total'] >= 1)]

df_2 = df_2.rename(columns={'gene_name': 'GN'})

df = df_2.merge(df_1, how='left', on=['GN'])
df_select = df.loc[df['verify'] == 1]


mg = mygene.MyGeneInfo()
df_2_gene_list = df_2.loc[:, 'gene_id'].values.tolist()
df_2_ezID_list = []
# parameter: ['ensembl.gene', 'symbol', 'entrezgene', 'uniprot']
for s in mg.querymany(df_2_gene_list, scopes='ensembl.gene', fields='entrezgene, symbol', species='human'):
    if 'entrezgene' in s:
        df_2_ezID_list.append({s['query']: s['entrezgene'], 'symbol': s['symbol']})

df_1_ezID_list = []
# df_1_ezID_dict = {}
df_1_gene_list = df_1.loc[:, 'GN'].values.tolist()
df_1_accessionID = df_1.loc[:, 'Accession'].values.tolist()
for q in mg.querymany(df_1_accessionID, scopes='uniprot', fields='entrezgene, symbol', species='human'):
    if 'entrezgene' in q:
        df_1_ezID_list.append({q['query']: q['entrezgene'], 'symbol': q['symbol']})
        # df_1_ezID_dict[q['query']] = q['entrezgene']
pprint(df_1_ezID_list)
print('='*100)
pprint(df_2_ezID_list)

with open(os.path.join('.', 'pickle_storage/total_pathway_info_list.pickle'), 'rb') as r:
    info_data_list = pickle.load(r)
print(info_data_list)

p1 = []
involved = []
for info in info_data_list:
    for i1 in df_1_ezID_list:
        if list(i1.values())[0] in info['involvedGeneID']:
            for i2 in df_2_ezID_list:
                if list(i2.values())[0] in info['involvedGeneID']:

                    temp = {'prot': f'{list(i1.keys())[0]}/{list(i1.values())[0]}',
                            'rnaseq': f'{list(i2.keys())[0]}/{list(i2.values())[0]}',
                            'metapathway': info['metaPathway'],
                            'meta_id': info['metaID'],
                            'subpathway': info['subPathway'],
                            'subpathway_id': info['subPathwayID'],
                            'prot_sym': i1['symbol'],
                            'prot_level': df_1['nor Heavy / Light_pho-total'].loc[df_1['Accession'] == list(i1.keys())[0]].values.tolist()[0],
                            'rnaseq_sym': i2['symbol'],
                            'rna_level': df_2['TPM變化量'].loc[df_2['gene_id'] == list(i2.keys())[0]].values.tolist()[0],
                            }

                    if temp not in involved:
                        involved.append(temp)
                    if info not in p1:
                        p1.append(info)

involved_df = pd.DataFrame(involved)
involved_df.to_excel(os.path.join('.', 'exported_data', 'rnaseq_prot_merge.xlsx'), index=False)

