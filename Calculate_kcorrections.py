
# Trying to get k-corrections for LCBG in cluster data. 
# Testing with Molino et al. catalog 2017 for cluster 
# MACS2129.
#
#Editor	Date			Change
#
# LH		10/12/2018		First write

import numpy as np
import pandas as pd
import kcorrect 
import kcorrect.utils as ut

infile='/home/lrhunt/Documents/LCBG_CLUSTER_GROUP/hlsp_clash_hst_ir_macs2129_cat-molino.txt'

dataset=pd.read_table(infile,delim_whitespace=True,header=133)

# Example of how to select rows in pandas

cluster_member=dataset[(dataset.zb_1<dataset.clusterz+0.2)&(dataset.zb_1>dataset.clusterz-0.2)]

# Example of how to select columns in pandas

cluster_member_photometry=cluster_member[['CLASHID','F435W_ACS_MASS','dF435W_ACS_MASS','F475W_ACS_MASS','dF475W_ACS_MASS','F606W_ACS_MASS','dF606W_ACS_MASS','F625W_ACS_MASS','dF625W_ACS_MASS','F775W_ACS_MASS','dF775W_ACS_MASS','F814W_ACS_MASS','dF814W_ACS_MASS','F850LP_ACS_MASS','dF850LP_ACS_MASS','F105W_WFC3_MASS','dF105W_WFC3_MASS','F110W_WFC3_MASS','dF110W_WFC3_MASS','F125W_WFC3_MASS','dF125W_WFC3_MASS','F140W_WFC3_MASS','dF140W_WFC3_MASS','F160W_WFC3_MASS','dF160W_WFC3_MASS','F225W_WFC3_MASS','dF225W_WFC3_MASS','F275W_WFC3_MASS','dF275W_WFC3_MASS','F336W_WFC3_MASS','dF336W_WFC3_MASS','F390W_WFC3_MASS','dF390W_WFC3_MASS','zb_1']]

cluster_member_maggies=cluster_member_photometry.copy(deep=True)

cluster_member_maggies=cluster_member_maggies[(cluster_member_maggies.F435W_ACS_MASS<80)&(cluster_member_maggies.F475W_ACS_MASS<80)&(cluster_member_maggies.F606W_ACS_MASS<80)&(cluster_member_maggies.F625W_ACS_MASS<80)&(cluster_member_maggies.F775W_ACS_MASS<80)&(cluster_member_maggies.F814W_ACS_MASS<80)&(cluster_member_maggies.F850LP_ACS_MASS<80)&(cluster_member_maggies.F105W_WFC3_MASS<80)&(cluster_member_maggies.F110W_WFC3_MASS<80)&(cluster_member_maggies.F125W_WFC3_MASS<80)&(cluster_member_maggies.F140W_WFC3_MASS<80)&(cluster_member_maggies.F160W_WFC3_MASS<80)&(cluster_member_maggies.F225W_WFC3_MASS<80)&(cluster_member_maggies.F275W_WFC3_MASS<80)&(cluster_member_maggies.F336W_WFC3_MASS<80)&(cluster_member_maggies.F390W_WFC3_MASS<80)]

cluster_member_maggies.iloc[:,np.arange(1,33,2)]=ut.mag2maggies(cluster_member_maggies.iloc[:,np.arange(1,33,2)])

for var in np.arange(2,33,2):
	cluster_member_maggies.iloc[:,var]=ut.invariance(cluster_member_maggies.iloc[:,var-1],cluster_member_maggies.iloc[:,var])

cluster_member_maggies["c1"]=np.nan
cluster_member_maggies["c2"]=np.nan
cluster_member_maggies["c3"]=np.nan
cluster_member_maggies["c4"]=np.nan
cluster_member_maggies["c5"]=np.nan
cluster_member_maggies["c6"]=np.nan


cluster_member_maggies_from_kcorr=cluster_member_maggies.iloc[:,np.arange(1,31,2)]
cluster_member_maggies_from_kcorr.insert(loc=0,column='redshift',value=np.nan)

cmm_ind=cluster_member_maggies.index.values

# May be able to create a new filter list based on the column names
# that are already loadable from the text file. If that is possible,
# don't have to make a bunch of different templates. May need to
# write a new file.

kcorrect.load_templates()
kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/LCBG_CLUSTER_FLITS.dat')

for i in cmm_ind:
	#cluster_member_maggies.loc[i,'c1':'c6']=kcorrect.fit_nonneg(np.array(cluster_member_maggies.loc[i,'zb_1'],dtype=float),np.array(cluster_member_maggies.loc[i, np.array(list(cluster_member_maggies))[np.arange(1,31,2)]],dtype=float),np.array(cluster_member_maggies.loc[i,np.array(list(cluster_member_maggies))[np.arange(1,31,2)]],dtype=float))

#Need to add column for redshift in cluster_member_maggies_from_kcorr

for i in cmm_ind:
	cluster_member_maggies_from_kcorr.loc[i]=kcorrect.reconstruct_maggies(cluster_member_maggies.loc[i,'c1':'c6'])

columns=list(cluster_member_maggies_from_kcorr)

new_columns={}

for column in columns:
    new_columns.update({column:column+'_s'})

cluster_member_maggies_from_kcorr=cluster_member_maggies_from_kcorr.rename(columns=new_columns)

