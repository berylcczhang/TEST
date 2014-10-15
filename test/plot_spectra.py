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
modname='SLUG_SPECTRA'
allmodel=['_1','_3','_4','_2']
alllab=['Default','A$_V$ det.','A$_V$ prob.','z=3']

#start a bunch of display options
ncol=1
nraw=4

# Plot integrated quantities 
plt.figure(1, figsize=(15,24))
plt.subplots_adjust(hspace=0.2, wspace=0.3, bottom=0.05, top=0.95, 
                    left=0.10, right=0.98)

#loop over models 
nmodel=len(allmodel)

medians=[]
wavmed=[]

for mm in range(nmodel):
    
    #load output 
    cluster=read_cluster(modname+allmodel[mm])
    integrated=read_integrated(modname+allmodel[mm])
    

    #define spect and wave 
    try:
        mywave=integrated.wl_ex
        myspec=integrated.spec_ex
    except AttributeError:
        mywave=integrated.wl
        myspec=integrated.spec

    nwave=len(mywave)    

    #reconstruct stats for the spectra 
    spectra_int_med=np.zeros(shape=(nwave))
    spectra_int_75=np.zeros(shape=(nwave))
    spectra_int_25=np.zeros(shape=(nwave))
    spectra_int_05=np.zeros(shape=(nwave))
    spectra_int_95=np.zeros(shape=(nwave))
    
    for ww in range(nwave):
        spectra_int_med[ww]=np.median(myspec[ww,0,:])
        spectra_int_05[ww]=np.percentile(myspec[ww,0,:],05)
        spectra_int_25[ww]=np.percentile(myspec[ww,0,:],25)
        spectra_int_75[ww]=np.percentile(myspec[ww,0,:],75)
        spectra_int_95[ww]=np.percentile(myspec[ww,0,:],95)
        
    #store medians
    medians.append(spectra_int_med)    
    wavmed.append(mywave)    

    #display     
    sbpn=plt.subplot(nraw,ncol,mm+1)
    plt.xlabel("Wavelength (A)")
    plt.ylabel("L (erg/s/A)")
    plt.title(alllab[mm])

    #5-95%
    plt.fill_between(mywave,spectra_int_05,spectra_int_95,facecolor='red', alpha=0.5)
    plt.plot(mywave,spectra_int_05,color='red', linewidth=3)
    plt.plot(mywave,spectra_int_95,color='red', linewidth=3)
    
    #25-75%
    plt.fill_between(mywave,spectra_int_25,spectra_int_75,facecolor='blue', alpha=0.6)
    plt.plot(mywave,spectra_int_25,color='blue', linewidth=3)
    plt.plot(mywave,spectra_int_75,color='blue', linewidth=3)

    #median
    plt.plot(mywave,spectra_int_med,color='black', linewidth=3 )

    plt.xlim(500,10000)
    plt.ylim(1e32,5e37)
    plt.loglog()



#check if Av is working     
#sbpn=plt.subplot(nraw,ncol,5)
#plt.xlabel("Wavelength (A)")
#plt.ylabel("L (erg/s/A)")
#plt.xlim(100,10000)
#plt.ylim(1e30,1e38)
#plt.loglog()
#for mm in range(nmodel):
#    plt.plot(wavmed[mm],medians[mm], linewidth=3)
    

# Save
plt.savefig('test/'+modname+'_f1.pdf')



