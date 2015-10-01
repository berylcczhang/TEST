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
modname='SLUG_IMFCHOICE'
allmodel=['_1','_2','_3']
alllab=['Kroupa','Chabrier','Salpeter']

#load output 
cluster= []

for mm in allmodel:
    cluster.append(read_cluster(modname+mm))


#number cluster 
ncls   = len(cluster[0].id)
ntim   = len(cluster[0].time)

#start a bunch of display options

#Plot cluster quantities 
plt.figure(1, figsize=(12,18))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.95)


#plot masses of cluster 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(4,3,loop+1)
    plt.hist(cluster[loop].max_star_mass,bins=20,range=(0,120))
    plt.title(alllab[loop],fontsize=15)
    plt.xlabel('Max Stellar Mass (M$_\odot$)')
    plt.ylabel('Number')
    plt.xlim(0,150)
    plt.ylim(0,250)
    plt.axvline(x=cluster[0].target_mass[0],color='red',linewidth=2.5,linestyle='--')
    loop=loop+1


#plot number of stars 
loop=0
for mm in allmodel:
    sbpn=plt.subplot(4,3,loop+4)
    plt.hist(np.log10(cluster[loop].num_star),bins=20,range=(2.5,3.5))
    plt.xlabel('Number stars')
    plt.ylabel('Number')
    plt.xlim(2.5,3.5)
    plt.ylim(0,400)
    loop=loop+1



#plot QH0 of cluster 
loop=0
indx=cluster[0].filter_names.index('QH0')
for mm in allmodel:
    sbpn=plt.subplot(4,3,loop+7)
    plt.hist(cluster[loop].phot[:,indx]/1e49,bins=20,range=(0,20))
    plt.xlabel('QH0'+' (10$^{49}$ '+cluster[loop].filter_units[indx]+')')
    plt.ylabel('Number')
    plt.xlim(0,20)
    plt.ylim(0,400)
    loop=loop+1


#plot FUV of cluster 
loop=0
indx=cluster[0].filter_names.index('Lbol')
for mm in allmodel:
    sbpn=plt.subplot(4,3,loop+10)
    plt.hist(cluster[loop].phot[:,indx]/1e5,bins=20,range=(0,20))
    plt.xlabel('Lbol'+' (10$^{5}$ '+cluster[loop].filter_units[indx]+')')
    plt.ylabel('Number')
    plt.xlim(0,25)
    plt.ylim(0,400)
    loop=loop+1

# Save
plt.savefig('test/'+modname+'_f1.pdf')


