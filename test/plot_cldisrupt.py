#
#This script displays the output of the test problem: example
#See section 10 of the slug manual 
#

import numpy as np
import matplotlib.pyplot as plt
import os
import os.path as osp
import copy
from os import listdir 

try:
    from slugpy import *    # If slugpy is already in our path
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        import sys
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Get slug dir
if 'SLUG_DIR' in os.environ:
    slugdir = os.environ['SLUG_DIR']
else:
    slugdir = '.'


#set model name and output type
modname='SLUG_CLDISRUPT'
allmodel=['_1','_2']
alllab=['No Disr','With Disr.']

#load output 
cluster= []
integrated= []

for mm in allmodel:
    cluster.append(read_cluster(modname+mm))
    integrated.append(read_integrated(modname+mm))

#start a bunch of display options
ncol=2
nraw=4

ntime=len(integrated[0].time)
nmodel=len(allmodel)

# Plot integrated quantities 
plt.figure(1, figsize=(12,20))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.98)

#plot masses in cluster 

for mm in range(nmodel):
    
    #find median, quartiles of mass in clusters 
    median_cluster=np.zeros(shape=(ntime))
    first_cluster=np.zeros(shape=(ntime))
    third_cluster=np.zeros(shape=(ntime))
    median_galaxy=np.zeros(shape=(ntime))
    median_ncluster=np.zeros(shape=(ntime))
    median_ndistcluster=np.zeros(shape=(ntime))
    
    median_lbol=np.zeros(shape=(ntime))
    first_lbol=np.zeros(shape=(ntime))
    third_lbol=np.zeros(shape=(ntime))
    median_fuv=np.zeros(shape=(ntime))
    first_fuv=np.zeros(shape=(ntime))
    third_fuv=np.zeros(shape=(ntime))


    #grab UV and Lbol
    indx_fuv=integrated[mm].filter_names.index('GALEX_FUV')
    indx_lbol=integrated[mm].filter_names.index('Lbol')

    for tt in range(ntime):
        median_cluster[tt]=np.median(integrated[mm].cluster_mass[tt,:])
        first_cluster[tt]=np.percentile(integrated[mm].cluster_mass[tt,:],25)
        third_cluster[tt]=np.percentile(integrated[mm].cluster_mass[tt,:],75)
        median_galaxy[tt]=np.median(integrated[mm].actual_mass[tt,:])
        median_ncluster[tt]=np.median(integrated[mm].num_clusters[tt,:])
        median_ndistcluster[tt]=np.median(integrated[mm].num_dis_clusters[tt,:])
        median_lbol[tt]=np.median(integrated[mm].phot[indx_lbol,tt,:])
        first_lbol[tt]=np.percentile(integrated[mm].phot[indx_lbol,tt,:],25)
        third_lbol[tt]=np.percentile(integrated[mm].phot[indx_lbol,tt,:],75)
        median_fuv[tt]=np.median(integrated[mm].phot[indx_fuv,tt,:])
        first_fuv[tt]=np.percentile(integrated[mm].phot[indx_fuv,tt,:],25)
        third_fuv[tt]=np.percentile(integrated[mm].phot[indx_fuv,tt,:],75)
     
    #plot cluster/galaxy mass
    sbpn=plt.subplot(nraw,ncol,mm+1)
    plt.title(alllab[mm],fontsize=15)
  
    time=integrated[mm].time/1e6

    plt.fill_between(time,first_cluster,third_cluster,facecolor='red', alpha=0.5)
    plt.plot(time,first_cluster,color='red', linewidth=3)
    plt.plot(time,third_cluster,color='red', linewidth=3)
    plt.plot(time,median_cluster,color='black', linewidth=3)
    plt.plot(time,median_galaxy,color='blue', linewidth=3, linestyle='--')

    plt.xlabel('Time (Myr)')
    plt.ylabel('Mass (M$_\odot$)')
    plt.yscale('log')
    plt.ylim(100,2000)
    
    #plot cluster/galaxy number
    sbpn=plt.subplot(nraw,ncol,mm+3)
  
    plt.plot(time,median_ncluster,color='black', linewidth=3)
    plt.plot(time,median_ndistcluster,color='blue', linewidth=3, linestyle='--')

    plt.xlabel('Time (Myr)')
    plt.ylabel('Number')
    plt.ylim(-2,18)
 
    #plot Lbol
    sbpn=plt.subplot(nraw,ncol,mm+5)
  
    plt.fill_between(time,first_lbol,third_lbol,facecolor='red', alpha=0.5)
    plt.plot(time,first_lbol,color='red', linewidth=3)
    plt.plot(time,third_lbol,color='red', linewidth=3)
    plt.plot(time,median_lbol,color='black', linewidth=3)

    plt.xlabel('Time (Myr)')
    plt.ylabel('Lbol'+' ('+integrated[mm].filter_units[indx_lbol]+')')
    plt.yscale('log')
    plt.ylim(1e5,1e6)

    #plot UV
    sbpn=plt.subplot(nraw,ncol,mm+7)
  
    plt.fill_between(time,first_fuv,third_fuv,facecolor='red', alpha=0.5)
    plt.plot(time,first_fuv,color='red', linewidth=3)
    plt.plot(time,third_fuv,color='red', linewidth=3)
    plt.plot(time,median_fuv,color='black', linewidth=3)

    plt.xlabel('Time (Myr)')
    plt.ylabel('FUV'+' ('+integrated[mm].filter_units[indx_fuv]+')')
    plt.yscale('log')
    plt.ylim(1e23,2e24)
    
  
# Save
plt.savefig('test/'+modname+'_f1.pdf')


