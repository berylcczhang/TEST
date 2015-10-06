import numpy as np
import matplotlib.pyplot as plt
import os
import os.path as osp
import copy
from subprocess import call
from scipy import stats
from scipy.special import erfinv
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

# IMFs to use
imfs = ['chabrier', 'kroupa', 'wk06']

# Get slug dir
if 'SLUG_DIR' in os.environ:
    slugdir = os.environ['SLUG_DIR']
else:
    slugdir = '.'

# Loop over imfs
fnames=[]
for imf in imfs:

    # Print status
    print("Running with IMF "+imf+"...")

    # Create param file
    fpin = open(osp.join(slugdir, 'param',
                         'example_cluster.param'), 'r')
    fpout = open(osp.join(slugdir, 'param', 
                          'cluster_imf_compare.param'), 'w')
    for line in fpin:
        linesplit = line.split()
        if len(linesplit) == 0:
            fpout.write(line)
        elif linesplit[0] == 'output_mode':
            fpout.write('output_mode    fits')
        elif linesplit[0] == 'model_name':
            fnames.append('IMF_COMPARE_'+imf)
            fpout.write('model_name     '+fnames[-1]+'\n')
        elif linesplit[0] == 'imf':
            fpout.write('imf            lib/imf/'+imf+'.imf')
        elif linesplit[0] == 'phot_mode':
            fpout.write('phot_mode      L_lambda')
        else:
            fpout.write(line)
    fpin.close()
    fpout.close()

    # Run simulation
    call(['python ' + osp.join('bin', 'slug.py') + ' ' + 
          osp.join('param', 'cluster_imf_compare.param')],
         shell=True)

# Read runs
tcomp = [1e6, 2e6, 4e6, 1e7]
data = []
for f in fnames:
    data.append(read_cluster('output/'+f))

# Get density map at each wavelength at each time
wl = data[0].wl
specden = np.zeros((len(fnames), len(tcomp), 1000, 1221))
for i, d in enumerate(data):
    for j, t in enumerate(tcomp):
        spec = d.spec[np.where(d.time == t)[0],:]
        for k in range(1221):
            idx = np.where(spec[:,k]>0.0)[0]
            if len(idx) == 0:
                continue
            spec1 = np.log10(spec[idx,k])
            xvec = np.linspace(np.amin(spec1), np.amax(spec1), 20)
            tmp = stats.gaussian_kde(spec1)(xvec)
            specden[i,j,idx,k] = np.interp(spec1, xvec, tmp)
            specden[i,j,idx,k] \
                = specden[i,j,idx,k]/np.amax(specden[i,j,idx,k])

# Get mean and percentiles for spectra
meanspec = np.zeros((len(fnames), len(tcomp), 1221))
percspec = np.zeros((len(fnames), len(tcomp), 3, 1221))
for i, d in enumerate(data):
    for j, t in enumerate(tcomp):
        spec = d.spec[np.where(d.time == t),:].reshape(1000, 1221)
        meanspec[i,j,:] = np.mean(spec, axis=0)
        percspec[i,j,:,:] = np.percentile(spec, (10, 50, 90), axis=0)

# Get mean and percentiles for photometry; also get errors on percentiles
filter_wl_eff = data[0].filter_wl_eff
filter_wl = data[0].filter_wl
filter_response = data[0].filter_response
QH0idx = data[0].filter_names.index('QH0')
meanphot = np.zeros((len(fnames), len(tcomp), len(data[0].filter_names)))
percphot = np.zeros((len(fnames), len(tcomp), 3, len(data[0].filter_names)))
names=['$Q(\mathrm{H}^0)$', '$L_{\mathrm{bol}}$', 'FUV', 'F225',
       'F336', 'F555', 'F814']
for i, d in enumerate(data):
    for j, t in enumerate(tcomp):
        phot = d.phot[np.where(d.time == t),:]. \
               reshape(1000, len(data[0].filter_names))
        meanphot[i,j,:] = np.mean(phot, axis=0)
        percphot[i,j,:,:] = np.percentile(phot, (10, 50, 90), axis=0)
        phot.sort(axis=0)
        ndat = phot.shape[0]
        rankerr = np.sqrt(2*0.9*0.1*ndat)*erfinv(0.68)
        rankerrmed = np.sqrt(2*0.5**2*ndat)*erfinv(0.68)
        print "Time {:f} Myr,".format(t/1e6) + " IMF " + imfs[i] + \
            ": 10th percentile     50th percentile    90th percentile"
        for k in range(phot.shape[1]):
            print "{:10s}: {:f} -{:f} +{:f} dex    {:f} -{:f} +{:f} dex  {:f} -{:f} +{:f} dex".\
                format(names[k], np.log10(phot[0.1*ndat,k]), 
                       np.log10(phot[0.1*ndat,k]/phot[0.1*ndat-np.ceil(rankerr),k]),
                       np.log10(phot[0.1*ndat+np.ceil(rankerr),k]/phot[0.1*ndat,k]),
                       np.log10(phot[0.5*ndat,k]), 
                       np.log10(phot[0.5*ndat,k]/phot[0.5*ndat-np.ceil(rankerrmed),k]),
                       np.log10(phot[0.5*ndat+np.ceil(rankerrmed),k]/phot[0.5*ndat,k]),
                       np.log10(phot[0.9*ndat,k]), 
                       np.log10(phot[0.9*ndat,k]/phot[0.9*ndat-np.ceil(rankerr),k]),
                       np.log10(phot[0.9*ndat+np.ceil(rankerr),k]/phot[0.9*ndat,k]))

# Get cumulative and differential distributions for photometry
nbin = 20
photsort = np.zeros((len(fnames), len(tcomp), 1000, len(data[0].filter_names)))
photbin = np.zeros((len(fnames), len(tcomp), nbin, len(data[0].filter_names)))
photbinedge = np.zeros((len(fnames), len(tcomp), nbin+1, len(data[0].filter_names)))
for i, d in enumerate(data):
    for j, t in enumerate(tcomp):
        phot = np.copy(d.phot[np.where(d.time == t),:]. \
                       reshape(1000, len(data[0].filter_names)))
        photsort[i,j,:,:] = np.sort(phot, axis=0)
        units = copy.deepcopy(data[0].filter_units)
        photometry_convert('AB', photsort[i,j,:,:], units,
                           data[0].filter_wl_eff, filter_last=True)
        for k in range(len(data[0].filter_names)):
            photbin[i,j,:,k], photbinedge[i,j,:,k] \
                = np.histogram(np.log10(phot[:,k]), bins=nbin)


# Plot spectra
fig = plt.figure(1, figsize=(10,8))
ptskip = 10
for i, t in enumerate(tcomp):

    # Create panel
    ax=plt.subplot(4,3,3*i+1)

    # Chabrier
    chabrier_mean,=plt.loglog(wl, meanspec[0,i,:], 'k', lw=2)
    spec = data[0].spec[np.where(data[0].time == t)[0],:]
    plt.scatter(np.tile(wl,(1000,))[::ptskip], np.ravel(spec)[::ptskip],
                marker='.', linewidths=0, 
                c=np.ravel(specden[0,i,:,:])[::ptskip],
                vmin=0, vmax=1, rasterized=True)

    # Axis labels
    if i != len(tcomp)-1:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.ylabel(r'$L_\lambda$ [erg s$^{-1}$ $\AA^{-1}$]')

    # Limits
    plt.ylim([2e33,2e37])
    plt.xlim([1e2, 2e4])

    # Title
    if i == 0:
        plt.title('Chabrier')

    # Create panel
    ax=plt.subplot(4,3,3*i+2)

    #  Kroupa
    kroupa_mean,=plt.loglog(wl, meanspec[1,i,:], 'k', lw=2)
    spec = data[1].spec[np.where(data[1].time == t)[0],:]
    plt.scatter(np.tile(wl,(1000,))[::ptskip], np.ravel(spec)[::ptskip],
                marker='.', linewidths=0, 
                c=np.ravel(specden[1,i,:,:])[::ptskip],
                vmin=0, vmax=1, rasterized=True)

    # Axis labels
    if i != len(tcomp)-1:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.setp(ax.get_yticklabels(), visible=False)

    # Limits
    plt.ylim([2e33,2e37])
    plt.xlim([1e2, 2e4])

    # Title
    if i == 0:
        plt.title('Kroupa')

    # Create panel
    ax=plt.subplot(4,3,3*i+3)

    #  WK06
    wk06_mean,=plt.loglog(wl, meanspec[2,i,:], 'k', lw=2)
    spec = data[2].spec[np.where(data[2].time == t)[0],:]
    plt.scatter(np.tile(wl,(1000,))[::ptskip], np.ravel(spec)[::ptskip],
                marker='.', linewidths=0, 
                c=np.ravel(specden[2,i,:,:])[::ptskip],
                vmin=0, vmax=1, rasterized=True)

    # Axis labels
    if i != len(tcomp)-1:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.setp(ax.get_yticklabels(), visible=False)

    # Limits
    plt.ylim([2e33,2e37])
    plt.xlim([1e2, 2e4])

    # Title
    if i == 0:
        plt.title('Weidner-Kroupa')

    # Labels
    plt.text(2.5e3, 4e36, 't = {:d} Myr'.format(int(tcomp[i]/1e6)))

# Adjust subplot spacing
plt.subplots_adjust(hspace=0, wspace=0, bottom=0.1, top=0.95, 
                    left=0.12, right=0.88)

# Add colorbar
axcbar = fig.add_axes((0.88, 0.1, 0.03, 0.85))
cbar = plt.colorbar(cax=axcbar)
cbar.set_label('Probability density')

# Save
plt.savefig('imfvary1.pdf')

# Plot photometry
offset=1.075
plt.figure(2, figsize=(7,8))

for i, t in enumerate(tcomp):

    # Create panel
    ax=plt.subplot(4,1,i+1)

    # Chabrier, GALEX and HST filters
    chabrier_mean,=plt.loglog(filter_wl_eff[2:], meanphot[0,i,2:], 'go')
    chabrier_med,=plt.loglog(filter_wl_eff[2:], percphot[0,i,1,2:], 'o',
                             mfc='w', mec='g', mew=2)
    chabrier_spec,=plt.loglog(wl, meanspec[0,i,:], 'g', lw=1)
    plt.errorbar(filter_wl_eff[2:], percphot[0,i,1,2:],
                 yerr=np.abs(percphot[0,i,::2,2:]-percphot[0,i,1,2:]),
                 fmt='none', ecolor='g', elinewidth=2)

    # Kroupa, GALEX and HST filters
    kroupa_mean,=plt.loglog(np.array(filter_wl_eff[2:])*offset, 
                            meanphot[1,i,2:], 'bo')
    kroupa_med,=plt.loglog(np.array(filter_wl_eff[2:])*offset, 
                           percphot[1,i,1,2:], 'o',
                           mfc='w', mec='b', mew=2)
    kroupa_spec,=plt.loglog(wl, meanspec[1,i,:], 'b', lw=1)
    plt.errorbar(np.array(filter_wl_eff[2:])*offset, percphot[1,i,1,2:],
                 yerr=np.abs(percphot[1,i,::2,2:]-percphot[1,i,1,2:]),
                 fmt='none', ecolor='b', elinewidth=2)

    # WK06, GALEX and HST filters
    wk06_mean,=plt.loglog(np.array(filter_wl_eff[2:])/offset, 
                            meanphot[2,i,2:], 'ro')
    wk06_med,=plt.loglog(np.array(filter_wl_eff[2:])/offset, 
                           percphot[2,i,1,2:], 'o',
                           mfc='w', mec='r', mew=2)
    wk06_spec,=plt.loglog(wl, meanspec[2,i,:], 'r', lw=1)
    plt.errorbar(np.array(filter_wl_eff[2:])/offset, percphot[2,i,1,2:],
                 yerr=np.abs(percphot[2,i,::2,2:]-percphot[2,i,1,2:]),
                 fmt='none', ecolor='r', elinewidth=2)

    # Legend
    if i==len(tcomp)-1:
        plt.legend([chabrier_spec, kroupa_spec, wk06_spec],
                   ['Chabrier', 'Kroupa', 'WK'],
                   loc='lower left')

    # Filter responses
    if i==0:
        for j in range(len(filter_wl)):
            if filter_wl[j] is None:
                continue
            plt.plot(filter_wl[j], 
                     filter_response[j]*1.0e35
                     /np.amax(filter_response[j]), 'k')
            plt.text(filter_wl_eff[j], 2e32, names[j],
                     horizontalalignment='center',
                     bbox=dict(facecolor='w', lw=0))

    # Limits
    plt.ylim([1e32,4e36])
    plt.xlim([5e2, 1e4])

    # Ionization threshold
    if i==0:
        plt.plot([912, 912], [1e32, 1e35], 'k--')

    # Labels
    plt.ylabel(r'$L_\lambda$ [erg s$^{-1}$ $\AA^{-1}$]')
    if i != len(tcomp)-1:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.text(4e3, 1e36, 't = {:d} Myr'.format(int(tcomp[i]/1e6)))

    # Ionizing luminosity
    axt = ax.twinx()
    plt.loglog([700], [meanphot[0,i,QH0idx]], 'gs')
    plt.plot([700], [percphot[0,i,1,QH0idx]], 's',
             mfc='w', mec='g', mew=2)
    plt.errorbar([700], [percphot[0,i,1,QH0idx]],
                 yerr=np.abs(percphot[0,i,::2,QH0idx].reshape(2,1)
                             -[percphot[0,i,1,QH0idx]]),
                 fmt='none', ecolor='g', elinewidth=2)
    plt.loglog([700*offset], [meanphot[1,i,QH0idx]], 'bs')
    plt.plot([700*offset], [percphot[1,i,1,QH0idx]], 's',
             mfc='w', mec='b', mew=2)
    plt.errorbar([700*offset], [percphot[1,i,1,QH0idx]],
                 yerr=np.abs(percphot[1,i,::2,QH0idx].reshape(2,1)
                             -[percphot[1,i,1,QH0idx]]),
                 fmt='none', ecolor='b', elinewidth=2)
    plt.loglog([700/offset], [meanphot[2,i,QH0idx]], 'rs')
    plt.plot([700/offset], [percphot[2,i,1,QH0idx]], 's',
             mfc='w', mec='r', mew=2)
    plt.errorbar([700/offset], [percphot[2,i,1,QH0idx]],
                 yerr=np.abs(percphot[2,i,::2,QH0idx].reshape(2,1)
                             -[percphot[2,i,1,QH0idx]]),
                 fmt='none', ecolor='r', elinewidth=2)
    plt.ylim([1e48,2e50])
    plt.xlim([5e2, 1e4])
    plt.ylabel(r'$Q(\mathrm{H}^0)$ [phot s$^{-1}$]')
    if i==0:
        plt.text(700, np.exp(np.log(1e48) + 
                             np.log(2e50/1e48) *
                             np.log(2e32/1e32) / np.log(4e36/1e32)),
                 names[0], horizontalalignment='center',
                 bbox=dict(facecolor='w', lw=0))


# Adjust subplot spacing
plt.subplots_adjust(hspace=0, wspace=0, bottom=0.1, top=0.95, 
                    left=0.15, right=0.85)

# Save
plt.savefig('imfvary2.pdf')

# Set of filters and times to plot
filt = ['QH0', 'GALEX_FUV', 'WFC3_UVIS_F225W', 'WFC3_UVIS_F336W']
labels = [r'$Q(\mathrm{H}^0)$', 'FUV', 'F225', 'F336']
linestyle = ['', '--', '-.']
tfilt = [1e6, 4e6, 1e7]

# Plot cumulative distributions
plt.figure(3, figsize=(7,8))
for i, f in enumerate(filt):

    idx = data[0].filter_names.index(f)

    # Create panel
    ax=plt.subplot(len(filt),1,i+1)

    # Loop over times
    for j, t in enumerate(tfilt):

        # Plot CDF for different IMFs
        tidx = tcomp.index(t)
        pch,=plt.plot(photsort[0,tidx,:,idx], np.arange(1000)/1000., 
                      'g'+linestyle[j])
        pk,=plt.plot(photsort[1,tidx,:,idx], np.arange(1000)/1000., 
                     'b'+linestyle[j])
        pwk06,=plt.plot(photsort[2,tidx,:,idx], np.arange(1000)/1000., 
                        'r'+linestyle[j])
        # Legend
        if i==0 and j==0:
            plt.legend([pch, pk, pwk06],
                       ['Chabrier', 'Kroupa', 'WK'],
                       loc='upper left')

    # Extra legend
    if i==1:
        pt=[]
        tlabel=[]
        for j, t in enumerate(tfilt):
            pttmp,=plt.plot([-100,-100], [100,110], 'k'+linestyle[j])
            pt.append(pttmp)
            tlabel.append("t = {:d} Myr".format(int(t/1e6)))
        plt.legend(pt, tlabel, loc='upper left')

    # Set up x axis
    if f != 'QH0':
        plt.xlim([-3,-9])
    else:
        plt.xscale('log')
        plt.xlim([1e44,1e51])

    # y limit
    if i == 0:
        plt.ylim([0,1])
    else:
        plt.ylim([0,0.99])

    # Labels
    plt.ylabel('CDF')
    if i == len(filt)-1:
        plt.xlabel("M$_{\mathrm{AB}}$")
    elif i != 0:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        plt.xlabel("$Q(\mathrm{H}^0)$")
    if i != 0:
        plt.text(-8, 0.2, labels[i], bbox=dict(facecolor='w'), 
                 horizontalalignment='center')
    else:
        plt.text(1e44*(1e51/1e44)**((-8.+3.)/(-9.+3.)), 0.2,
                 labels[i], horizontalalignment='center',
                 bbox=dict(facecolor='w'))


# Adjust subplot spacing
plt.subplots_adjust(hspace=0, wspace=0, bottom=0.1, top=0.9, 
                    left=0.1, right=0.95)

# Save
plt.savefig('imfvary3.pdf')
