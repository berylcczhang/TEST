# This script processes the data files produced by make_grid into the
# format that SLUG reads.

import os
import os.path as osp
from subprocess import call
import shutil
import numpy as np
from slugpy import read_cluster_phot
from slugpy import read_integrated_phot
from slugpy.cloudy import read_cluster_cloudylines
from slugpy.cloudy import read_integrated_cloudylines
from slugpy.int_tabulated import int_tabulated

# List of tracks and their metallicities relative to solar
track_list = [
    # MIST non-rotating
    'mist/vvcrit000/MIST_v1.0_feh_m4.00_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m3.50_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m3.00_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m2.50_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m2.00_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m1.75_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m1.50_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m1.25_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m1.00_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m0.75_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m0.50_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_m0.25_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_p0.25_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    'mist/vvcrit000/MIST_v1.0_feh_p0.50_afe_p0.0_vvcrit0.0_EEPS.fits.gz',
    # MIST rotating
    'mist/vvcrit040/MIST_v1.0_feh_m4.00_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m3.50_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m3.00_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m2.50_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m2.00_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m1.75_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m1.50_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m1.25_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m1.00_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m0.75_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m0.50_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_m0.25_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_p0.25_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    'mist/vvcrit040/MIST_v1.0_feh_p0.50_afe_p0.0_vvcrit0.4_EEPS.fits.gz',
    # Geneva 2013 non-rotating models
    'sb99/Z0140v00.txt', 'sb99/Z0020v00.txt', 
    # Padova with TP-AGB models
    'sb99/modp0004.dat', 'sb99/modp004.dat', 'sb99/modp008.dat',
    'sb99/modp020.dat', 'sb99/modp050.dat',
    # Geneva 2014 rotating models
    'Z0020v40.txt', 'Z0140v40.txt',
    # Geneva pre-2013, standard mass loss
    'sb99/modc001.dat', 'sb99/modc004.dat', 'sb99/modc008.dat',
    'sb99/modc020.dat', 'sb99/modc040.dat',
    # Geneva pre-2013, enhanced mass los
    'sb99/mode001.dat', 'sb99/mode004.dat', 'sb99/mode008.dat',
    'sb99/mode020.dat', 'sb99/mode040.dat',
    # Padova w/out TP-AGB
    'sb99/mods0004.dat', 'sb99/mods004.dat', 'sb99/mods008.dat',
    'sb99/mods020.dat', 'sb99/mods050.dat'
]

# List of ionization parameters
logU_list = [-3, -2.5, -2]

# Name of template cloudy file
cloudy_input = osp.join(os.environ['SLUG_DIR'],
                        'cloudy_slug', 'cloudy.in_grid_template')

# Name of the SLUG parameter files
slug_param = osp.join(os.environ['SLUG_DIR'], 'param',
                      'cluster_cts_10myr.param')
slug_param_gal = osp.join(os.environ['SLUG_DIR'],
                          'param', 'galaxy_cts_10myr.param')

# Read files
alldata = []
for track in track_list:

    # Read the input files to get the model name
    slug_param_tmp = '.'.join(slug_param.split('.')[:-1]) + \
                     '_' + '.'.join(
                         osp.basename(track).split('.')[:-1]) + \
                     '.param'
    fpin = open(slug_param, 'r')
    fpout = open(slug_param_tmp, 'w')
    for line in fpin:
        spl = line.split()
        if len(spl) > 0:
            if spl[0] == 'model_name':
                line = line[:-1] + '_' + '.'.join(
                    osp.basename(track).split('.')[:-1])
                modname = line.split()[-1]
            elif spl[0] == 'tracks':
                line = 'tracks   lib/tracks/' + track
        fpout.write(line)
    fpin.close()
    fpout.close()
    slug_param_tmp_gal = '.'.join(slug_param_gal.split('.')[:-1]) + \
                     '_' + '.'.join(
                         osp.basename(track).split('.')[:-1]) + \
                     '.param'
    fpin = open(slug_param_gal, 'r')
    fpout = open(slug_param_tmp_gal, 'w')
    for line in fpin:
        spl = line.split()
        if len(spl) > 0:
            if spl[0] == 'model_name':
                line = line[:-1] + '_' + '.'.join(
                    osp.basename(track).split('.')[:-1])
                modname_gal = line.split()[-1]
            elif spl[0] == 'tracks':
                line = 'tracks   lib/tracks/' + track
        fpout.write(line)
    fpin.close()
    fpout.close()

    # Loop over ionization parameters
    for logU in logU_list:

        # Read the integrated data
        print "Reading {:s}, log U = {:f}...".format(modname_gal, logU)
        basename_gal = modname_gal+'_{:f}'.format(logU)
        dataphot_gal = read_integrated_phot(osp.join(
            "nebular_grid", basename_gal))
        datalines_gal = read_integrated_cloudylines(osp.join(
            "nebular_grid", basename_gal))

        # Read the cluster data
        print "Reading {:s}, log U = {:f}...".format(modname, logU)
        basename = modname+'_{:f}'.format(logU)
        dataphot = read_cluster_phot(osp.join(
            "nebular_grid", basename))
        datalines = read_cluster_cloudylines(osp.join(
            "nebular_grid", basename))

        # Read hydrogen conditions file for the galaxy run, and extract
        # the mean temperature
        fname = osp.join('cloudy_tmp_'+basename_gal,
                         basename_gal+'_tr00000_ti00000.hcon')
        hcon = np.loadtxt(fname)
        r = hcon[:,0]
        te = hcon[:,1]
        nhplus = hcon[:,2]*hcon[:,5]
        tgal = int_tabulated(r, 4.0*np.pi*r**2*nhplus**2*te) / \
               int_tabulated(r, 4.0*np.pi*r**2*nhplus**2)

        # Get temperature versus time for the cluster runs
        ntime = len(dataphot.time)
        tcl = np.zeros(ntime)
        for i in range(ntime):
            fname = osp.join('cloudy_tmp_'+basename,
                             basename+'_n{:09d}.hcon'.format(i))
            hcon = np.loadtxt(fname)
            r = hcon[:,0]
            te = hcon[:,1]
            nhplus = hcon[:,2]*hcon[:,5]
            tcl[i] = int_tabulated(r, 4.0*np.pi*r**2*nhplus**2*te) / \
                     int_tabulated(r, 4.0*np.pi*r**2*nhplus**2)

        # Read the line array to get the exact wavelengths
        linearr_name = osp.join('cloudy_tmp_'+basename,
                                basename+'_n000000000.linearr')
        fp = open(linearr_name, 'r')
        allwl = {}
        for i, line in enumerate(fp):
            spl = line.split('\t')
            try:
                lab = spl[1]
                wl = float(spl[0])
                allwl[lab] = wl
            except:
                pass

        # Construct a label in cloudy format to extract from the line
        # list
        linelabel = []
        for lab, wl in zip(datalines.cloudy_linelist, 
                           datalines.cloudy_linewl):
            if wl < 1e2:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:6.4f}A'.format(wl))
            elif wl < 1e3:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:6.3f}A'.format(wl))
            elif wl < 1e4:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:6.2f}A'.format(wl))
            elif wl < 1e5:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:7.5f}m'.format(wl/1e4))
            elif wl < 1e6:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:7.4f}m'.format(wl/1e4))
            else:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '      {:7.3f}m'.format(wl/1e4))
        linelabel = np.array(linelabel)

        # Identify entries that are not lines in the line array, and
        # throw these out
        idxkeep = []
        for i, lab in enumerate(linelabel):
            if lab in allwl:
                idxkeep.append(i)
        linelabel = linelabel[idxkeep]
        linelum = datalines.cloudy_linelum[:,idxkeep]
        linelum_gal = np.squeeze(datalines_gal.cloudy_linelum[idxkeep])[:,-1]

        # Identify H lines and throw these out
        idxkeep = []
        for i, lab in enumerate(linelabel):
            if lab[:2] != 'H ':
                idxkeep.append(i)
        linelabel = linelabel[idxkeep]
        linelum = linelum[:,idxkeep]
        linelum_gal = linelum_gal[idxkeep]

        # Normalize all line luminosities to ionizing luminosity
        lumnorm = np.transpose(linelum)/dataphot.phot[:,0]
        lumnorm_gal = linelum_gal/dataphot_gal.phot[0,-1,0]

        # Identify lines that are always below 10^-20 / ionizing
        # photon in luminosity, and throw these out
        idxkeep = np.logical_or(lumnorm_gal > 1e-20, 
                                np.sum(lumnorm > 1e-20, axis=1) > 0)
        linelabel = linelabel[idxkeep]
        lumnorm = lumnorm[idxkeep,:]
        lumnorm_gal = lumnorm_gal[idxkeep]

        # Get wavelengths for the remaining lines
        outwl=[]
        for lab in linelabel:
            outwl.append(allwl[lab])
        outwl = np.array(outwl)

        # Sort in order of decreasing luminosity at time 0
        idx = np.argsort(lumnorm[:,0])
        lumnorm = lumnorm[idx[::-1],:]
        lumnorm_gal = lumnorm_gal[idx[::-1]]
        outwl = outwl[idx[::-1]]
        linelabel = linelabel[idx[::-1]]

        # Append to master list
        outdict = { 
            'track' : osp.basename(track),
            'time' : dataphot.time,
            'logU' : logU,
            'wl' : outwl,
            'lum' : lumnorm,
            'lum_gal' : lumnorm_gal,
            'label' : linelabel,
            'tgal' : tgal,
            'tcl' : tcl 
        }
        alldata.append(outdict)

        #print "logU = {:f}, model = {:s}, Ha/QH0 = {:e}, T0 = {:f}, nline = {:d}".format(
        #    logU, modname, 
        #    datalines.cloudy_linelum[0,53]/dataphot.phot[0,0],
        #    tcl[0], len(linelabel))

# Write out the metal lines file header
fp = open(osp.join(os.environ['SLUG_DIR'], 'lib', 'atomic', 'cloudy_tab.txt'), 'w')
fp.write('# Metal lines and mean temperatures obtained by using SLUG to generate a\n')
fp.write('# fully-sampled spectrum for a Chabrier (2005) IMF at a range of ages,\n')
fp.write('# then using that spectrum in cloudy with constant density, background\n')
fp.write('# cosmic rays, expanding spehre geometry, and HII region abundances,\n')
fp.write('# scaled by  metallicity relative to Solar (which is taken to be 0.014\n')
fp.write('# for Geneva 2013 models, 0.020 for other models). The lines listed\n')
fp.write('# below are those for which the luminosity per ionizing photon\n')
fp.write('# injected is > 10^-20 erg/s/photon\n')
fp.write('#\n')
fp.write('# The first line below gives:\n')
fp.write('# NTRACK NTIME NLOGU\n')
fp.write('# where NTRACK is the number of tracks, NTIME is the number of\n')
fp.write('# output times, and NLOGU is the number of LOGU values. The next line\n')
fp.write('# is NTIME entries of the form\n')
fp.write('# TIME_1 TIME_2 TIME_3 ...\n')
fp.write('# where TIME_N is the nth output time. The line after that is\n')
fp.write('# LOGU_1 LOGU_2 LOGU_3 ...\n')
fp.write('# where LOGU_N is the Nth value of the ionization parameter. This is\n')
fp.write('# followed by a series of NTRACKS * NU blocks, each of which begins:\n')
fp.write('# TRACK LOGU NLINES\n')
fp.write('# where TRACK is the track name, LOGU is the log of the ionization\n')
fp.write('# parameter, and NLINES is the number of lines for this block. This is\n')
fp.write('# followed by one line of the form:\n')
fp.write('# TEMP_CTS   TEMP_1   TEMP_2   TEMP_3 ...\n')
fp.write('# where TEMP_CTS is the nH+^2-weighted mean temperature evaluated for\n')
fp.write('# the continuous star formation case, and TEMP_N is the temperature\n')
fp.write('# at each time for a single-age stellar population. This is followed\n')
fp.write('# by NLINES lines of the form\n')
fp.write('# SPEC     LAMBDA     (LUM_QH0)_CTS   (LUM/QH0)_1     (LUM/QH0)_2 ...\n')
fp.write('# where SPEC is the 11-character cloudy designation for the emitting\n')
fp.write('# species, LAMBDA is the line wavelength in Angstrom, (LUM/QH0)_CTS\n')
fp.write('# is the line luminosity per ionizing photon injected for the case.\n')
fp.write('# of continuous star formation, and (LUM/QH0)_N is the corresponding\n')
fp.write('# quantity at the Nth time for the single-age case.\n')
fp.write('\n')

# Write out the data header
fp.write('{:d}  {:d}  {:d}\n'.format(len(track_list), len(alldata[0]['time']),
                                     len(logU_list)))
for t in alldata[0]['time']:
    fp.write('{:e}   '.format(t))
fp.write('\n')
for logU in logU_list:
    fp.write('{:f}   '.format(logU))
fp.write('\n')

# Write out the data
for dat in alldata:
    fp.write('{:s}   {:f}   {:d}\n'.format(dat['track'], dat['logU'],
                                           len(dat['label'])))
    fp.write('{:e}'.format(dat['tgal']))
    for t in dat['tcl']:
        fp.write('   {:e}'.format(t))
    fp.write('\n')
    for label, wl, lum_gal, lum in \
        zip(dat['label'], dat['wl'], dat['lum_gal'], dat['lum']):
        fp.write('{:s}   {:e}   {:e}'.format(label, wl, lum_gal))
        for l in lum:
            fp.write('   {:e}'.format(l))
        fp.write('\n')

# Close the file
fp.close()
