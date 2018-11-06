
# coding: utf-8

# Below will be some code to calculate k-corrections for clusters in the Molino HST Cluster catalogs, following their example

# In[1]:


import numpy as np
import pandas as pd
import kcorrect 
import kcorrect.utils as ut
import matplotlib.pyplot as plt


# The next line reads the catalog into a pandas dataframe. Whitespace=True sets the delimiter (separator) as whitespace, and header=133 sets the line that contains header information. This line is the one that labels each column correctly in this document. The final line selects objects that have a photometric redshift within the cluster redshift range

# In[2]:


infile='/Users/lucashunt/projects/LCBGS_IN_CLUSTERS/hlsp_clash_hst_ir_macs2129_cat-molino.txt'
dataset=pd.read_table(infile,delim_whitespace=True,header=133)
cluster_member=dataset[(dataset.zb_1<dataset.clusterz+0.2)&(dataset.zb_1>dataset.clusterz-0.2)]


# This line selects only ID, photometry, and redshift columns

# In[3]:


cluster_member_photometry=cluster_member[['CLASHID','F435W_ACS_MASS','dF435W_ACS_MASS','F475W_ACS_MASS','dF475W_ACS_MASS','F606W_ACS_MASS','dF606W_ACS_MASS','F625W_ACS_MASS','dF625W_ACS_MASS','F775W_ACS_MASS','dF775W_ACS_MASS','F814W_ACS_MASS','dF814W_ACS_MASS','F850LP_ACS_MASS','dF850LP_ACS_MASS','F105W_WFC3_MASS','dF105W_WFC3_MASS','F110W_WFC3_MASS','dF110W_WFC3_MASS','F125W_WFC3_MASS','dF125W_WFC3_MASS','F140W_WFC3_MASS','dF140W_WFC3_MASS','F160W_WFC3_MASS','dF160W_WFC3_MASS','F225W_WFC3_MASS','dF225W_WFC3_MASS','F275W_WFC3_MASS','dF275W_WFC3_MASS','F336W_WFC3_MASS','dF336W_WFC3_MASS','F390W_WFC3_MASS','dF390W_WFC3_MASS','zb_1']]


# In[4]:


cluster_member_photometry=cluster_member_photometry[(cluster_member_photometry.F435W_ACS_MASS<80)&(cluster_member_photometry.F475W_ACS_MASS<80)&(cluster_member_photometry.F606W_ACS_MASS<80)&(cluster_member_photometry.F625W_ACS_MASS<80)&(cluster_member_photometry.F775W_ACS_MASS<80)&(cluster_member_photometry.F814W_ACS_MASS<80)&(cluster_member_photometry.F850LP_ACS_MASS<80)&(cluster_member_photometry.F105W_WFC3_MASS<80)&(cluster_member_photometry.F110W_WFC3_MASS<80)&(cluster_member_photometry.F125W_WFC3_MASS<80)&(cluster_member_photometry.F140W_WFC3_MASS<80)&(cluster_member_photometry.F160W_WFC3_MASS<80)&(cluster_member_photometry.F225W_WFC3_MASS<80)&(cluster_member_photometry.F275W_WFC3_MASS<80)&(cluster_member_photometry.F336W_WFC3_MASS<80)&(cluster_member_photometry.F390W_WFC3_MASS<80)]


# Below create a copy to generate maggies and select objects that actually have photometry

# In[5]:


cluster_member_maggies=cluster_member_photometry.copy(deep=True)


# In[6]:


cluster_member_maggies.iloc[:,np.arange(1,33,2)]=ut.mag2maggies(cluster_member_maggies.iloc[:,np.arange(1,33,2)])


# In[7]:


for var in np.arange(2,33,2):
	cluster_member_maggies.iloc[:,var]=ut.invariance(cluster_member_maggies.iloc[:,var-1],cluster_member_maggies.iloc[:,var])


# In[8]:


cluster_member_maggies["c1"]=np.nan
cluster_member_maggies["c2"]=np.nan
cluster_member_maggies["c3"]=np.nan
cluster_member_maggies["c4"]=np.nan
cluster_member_maggies["c5"]=np.nan
cluster_member_maggies["c6"]=np.nan


# In[9]:


cluster_member_maggies_from_kcorr=cluster_member_maggies.iloc[:,np.arange(1,31,2)]
cluster_member_maggies_from_kcorr.insert(loc=0,column='redshift',value=np.nan)


# In[10]:


cmm_ind=cluster_member_maggies.index.values


# In[11]:


kcorrect.load_templates()
kcorrect.load_filters('/Users/lucashunt/programs/kcorrect/data/templates/LCBG_CLUSTER_FLITS.dat')


# In[ ]:


for i in cmm_ind:
	cluster_member_maggies.loc[i,'c1':'c6']=kcorrect.fit_nonneg(np.array(cluster_member_maggies.loc[i,'zb_1'],dtype=float),np.array(cluster_member_maggies.loc[i, np.array(list(cluster_member_maggies))[np.arange(1,31,2)]],dtype=float),np.array(cluster_member_maggies.loc[i,np.array(list(cluster_member_maggies))[np.arange(1,31,2)]],dtype=float))


# In[ ]:


for i in cmm_ind:
	cluster_member_maggies_from_kcorr.loc[i]=kcorrect.reconstruct_maggies(cluster_member_maggies.loc[i,'c1':'c6'])


# In[ ]:


kcorrect.load_templates()
kcorrect.load_filters('/Users/lucashunt/programs/kcorrect/data/templates/bessell_ubv.dat')


# In[ ]:


bessel_ubv=pd.DataFrame({"U":[],"B":[],"V":[]})


# In[ ]:


print(len(cluster_member_maggies))
print(kcorrect.reconstruct_maggies(cluster_member_maggies[['c1','c2','c3','c4','c5','c6']].iloc[1]))


# In[ ]:


columns=list(cluster_member_maggies_from_kcorr)

new_columns={}

for column in columns:
    new_columns.update({column:column+'_s'})

cluster_member_maggies_from_kcorr=cluster_member_maggies_from_kcorr.rename(columns=new_columns)


# In[ ]:


print(cmm_ind)
print(-2.5*np.log10(cluster_member_maggies_from_kcorr.iloc[:,np.arange(1,len(new_columns))]))
print(cluster_member_photometry)


# In[ ]:


plt.plot(-2.5*np.log10(cluster_member_maggies_from_kcorr.iloc[:,3]),cluster_member_photometry.iloc[:,21],'.')

