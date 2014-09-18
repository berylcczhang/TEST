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
modname='SLUG_CLUSTER_EXAMPLE'

#load output 
cluster=read_cluster(modname)

#number cluster 
ncls   = len(cluster.id)
ntim   = len(cluster.time)

#start a bunch of display options

#Plot cluster quantities 
plt.figure(1, figsize=(8,15))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.95)

#grap last/first time step 
last=np.where(cluster.time == cluster.time.max())
first=np.where(cluster.time == cluster.time.min())

#plot masses of cluster 
sbpn=plt.subplot(4,2,1)
plt.hist(cluster.actual_mass[first[0]],bins=20,range=[480,520])
plt.xlim(480,520)
plt.xlabel('Actual Mass (M$_\odot$)')
plt.ylabel('Number')
plt.title('Time (Myr): '+str(cluster.time.min()/1e6),fontsize=20)

sbpn=plt.subplot(4,2,2)
plt.hist(cluster.actual_mass[last[0]],bins=20,range=[480,520])
plt.xlim(480,520)
plt.xlabel('Actual Mass (M$_\odot$)')
plt.ylabel('Number')
plt.title('Time (Myr): '+str(cluster.time.max()/1e6),fontsize=20)


#plot masses of most massive star alive 
sbpn=plt.subplot(4,2,3)
plt.hist(cluster.max_star_mass[first[0]])
plt.xlim(0,120)
plt.xlabel('Max Stellar Mass (M$_\odot$)')
plt.ylabel('Number')

sbpn=plt.subplot(4,2,4)
plt.hist(cluster.max_star_mass[last[0]])
plt.xlim(5,20)
plt.xlabel('Max Stellar Mass (M$_\odot$)')
plt.ylabel('Number')


#plot integrated photometry in selected bands
indx=cluster.filter_names.index('QH0')

sbpn=plt.subplot(4,2,5)
plt.hist(np.log10(cluster.phot[first[0],indx]))
plt.xlim(44,51)
plt.xlabel('log QH0'+' ('+cluster.filter_units[indx]+')')
plt.ylabel('Number')


sbpn=plt.subplot(4,2,6)
plt.hist(np.log10(cluster.phot[last[0],indx]))
plt.xlim(44,51)
plt.xlabel('log QH0'+' ('+cluster.filter_units[indx]+')')
plt.ylabel('Number')

#plot integrated photometry in selected bands
indx=cluster.filter_names.index('GALEX_FUV')

sbpn=plt.subplot(4,2,7)
plt.hist(np.log10(cluster.phot[first[0],indx]))
plt.xlim(21.5,24.5)
plt.xlabel('log GALEX_FUV'+' ('+cluster.filter_units[indx]+')')
plt.ylabel('Number')

sbpn=plt.subplot(4,2,8)
plt.hist(np.log10(cluster.phot[last[0],indx]))
plt.xlim(21.5,24.5)
plt.xlabel('log GALEX_FUV'+' ('+cluster.filter_units[indx]+')')
plt.ylabel('Number')

# Save
plt.savefig('test/'+modname+'_f1.pdf')


