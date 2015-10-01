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
modname='SLUG_CLFRACTION'
allmodel=['_1','_2','_3']
alllab=['fc=1.0','fc=0.5','fc=0.01']

#load output 
cluster= []
integrated= []

for mm in allmodel:
    cluster.append(read_cluster(modname+mm))
    integrated.append(read_integrated(modname+mm))

#start a bunch of display options
ncol=3
nraw=5

# Plot integrated quantities 
plt.figure(1, figsize=(12,20))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.98)


#plot masses of cluster 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(nraw,ncol,loop+1)
    plt.title(alllab[loop],fontsize=15)
    plt.hist(cluster[loop].max_star_mass,bins=20,range=(0,120))
    plt.xlabel('Max Stellar Mass in clusters (M$_\odot$)')
    plt.ylabel('Number')
    plt.xlim(0,150)
    plt.ylim(1,5000)
    plt.yscale('log')
    loop=loop+1

#plot number of stars in cluster 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(nraw,ncol,loop+4)
    plt.hist(cluster[loop].num_star,bins=20,range=(0,1000))
    plt.xlabel('Stars in clusters')
    plt.ylabel('Number')
    plt.xlim(0,1000)
    plt.ylim(1,1e4)
    plt.yscale('log')
    plt.text(400,1e3, 'N$_{total}$: '+str(sum(cluster[loop].num_star)))
    loop=loop+1

#plot field stars
loop=0
for mm in allmodel:
    sbpn=plt.subplot(nraw,ncol,loop+7)
    plt.hist(integrated[loop].num_fld_stars[0,:],bins=20,range=(1000,4000))
    plt.xlabel('Number Field Stars')
    plt.ylabel('Number')
    plt.xlim(1000,4000)
    plt.ylim(1,1e3)
    plt.text(2500,600, 'N$_{total}$: '+str(sum(integrated[loop].num_fld_stars[0,:])))
    plt.yscale('log')
    loop=loop+1

#plot target vs actual mass 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(nraw,ncol,loop+10)
    plt.hist(integrated[loop].actual_mass[0,:],bins=20,range=(0,3000))
    plt.xlabel('Actual Galaxy Mass (M$_\odot$)')
    plt.ylabel('Number')
    plt.xlim(0,3000)
    plt.ylim(1,1e3)
    plt.yscale('log')
    plt.axvline(x=integrated[0].target_mass[0,0],color='red',linewidth=2.5,linestyle='--')
    loop=loop+1



#plot photometry 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(nraw,ncol,loop+13)
    indx=integrated[loop].filter_names.index('QH0')
    plt.hist(integrated[loop].phot[indx,0,:]/1e49,bins=20,range=(0,35))
    plt.xlabel('QH0 (10$^{49}$ erg s$^{-1}$)')
    plt.ylabel('Number')
    plt.xlim(0,35)
    plt.ylim(1,1e2)
    plt.yscale('log')
    plt.axvline(x=np.median(integrated[loop].phot[indx,0,:]/1e49),color='red',linewidth=2.5,linestyle='--')
    loop=loop+1
      

# Save
plt.savefig('test/'+modname+'_f1.pdf')


