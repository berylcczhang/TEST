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
modname='SLUG_SFHSAMPLING'
allmodel=['_1','_2']
alllab=['fc=1.0','fc=0.0']

#load output 
cluster= []
integrated= []

for mm in allmodel:
    cluster.append(read_cluster(modname+mm))
    integrated.append(read_integrated(modname+mm))

#grab the num of time steps
ntime=len(integrated[0].time)
ntrial=len(integrated[0].actual_mass[0,:])
nmodel=len(allmodel)


#reconstruct the star formation histories 
sfharray=np.zeros(shape=(nmodel,ntrial,ntime))

for mm in range(nmodel):
    for tri in range(ntrial): 
        shift=np.roll(integrated[mm].actual_mass[:,tri],1)
        shift[0]=0
        orig=integrated[mm].actual_mass[:,tri]
        sfharray[mm,tri,:]=(orig-shift)/(integrated[mm].time[1]-integrated[mm].time[0])/1e-4
            
#start a bunch of display options
ncol=1
nraw=2

# Plot integrated quantities 
plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.98)

for mm in range(nmodel):
    sbpn=plt.subplot(nraw,ncol,mm+1)
    bp=plt.boxplot(sfharray[mm,:,:],showfliers=False,showmeans=True,whis=[5,95])
    plt.xlabel("Time Step")
    plt.ylabel("SFR (10$^{-4}$ M$_\odot$ yr$^{-1}$)")
    plt.title(alllab[mm])
    
    #some style
    for box in bp['boxes']:
        # change outline color
        box.set(linewidth=3)
    for whisker in bp['whiskers']:
        whisker.set(linewidth=3)
    for cap in bp['caps']:
        cap.set(linewidth=3)
    for median in bp['medians']:
        median.set(linewidth=3)

# Save
plt.savefig('test/'+modname+'_f1.pdf')


