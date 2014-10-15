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
modname='SLUG_GALAXY_EXAMPLE'

#load output 
cluster=read_cluster(modname)
integrated=read_integrated(modname)

#number time steps 
ntime = len(integrated.time)
nmod = len(integrated.actual_mass[0,:])

#start a bunch of display options

# Plot integrated quantities 
plt.figure(1, figsize=(8,12))
plt.subplots_adjust(hspace=0.2, wspace=0.25, bottom=0.05, top=0.95, 
                    left=0.10, right=0.98)



#mass in realization vs mass expected 
sbpn=plt.subplot(3,2,1)
for mm in range(nmod):
    plt.scatter(integrated.target_mass[:,mm]/1e6,integrated.actual_mass[:,mm]/1e6)
    plt.plot(np.arange(0,2,0.1),np.arange(0,2,0.1),color='red',linewidth=2.,linestyle='--')
    plt.xlim(0.00,0.23)
    plt.ylim(0.00,0.23)
    plt.xlabel('Target Mass ($10^6$ M$_\odot$)')
    plt.ylabel('Actual Mass ($10^6$ M$_\odot$)')

#plot integrated photometry in selected bands for 1 timestep 
select_filters=['QH0','Lbol','WFC3_UVIS_F225W','WFC3_UVIS_F555W','WFC3_UVIS_F814W']
subpan=2

for fil in integrated.filter_names:
    if fil in select_filters:
        sbpn=plt.subplot(3,2,subpan)
        indx=integrated.filter_names.index(fil)
        plt.scatter(integrated.actual_mass/1e6,integrated.phot[indx,:,:])
        plt.xlim(0.00,0.22)
        plt.xlabel('Actual Mass ($10^6$ M$_\odot$)')
        plt.ylabel(fil+' ('+integrated.filter_units[indx]+')')
        subpan=subpan+1
        

#for Ha, add exptected value given SFR 
sbpn=plt.subplot(3,2,2)
sfr=0.001/7.638e-54
plt.plot([-1,2],[sfr,sfr],color='red',linewidth=2.,linestyle='--')

#for Lbol, add exptected value given SFR 
sbpn=plt.subplot(3,2,4)
sfr=0.001/2.661e-44/3.846e33
plt.plot([-1,2],[sfr,sfr],color='red',linewidth=2.,linestyle='--')

# Save
plt.savefig('test/'+modname+'_f1.pdf')


