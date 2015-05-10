# This script generates the metal line grid that SLUG uses for its
# "quick and dirty" estimate of metal line luminosities.

import os
import os.path as osp
from subprocess import call
import shutil
from slugpy import read_cluster_phot
from slugpy.cloudy import read_cluster_cloudylines

# List of tracks and their metallicities relative to solar
track_list = ['Z0140v00.txt', 'Z0020v00.txt', 'modp0004.dat', 
              'modp004.dat', 'modp008.dat', 'modp020.dat',
              'modp050.dat' ]
zrel_list = [ 1.0, 1.0/7.0, 0.02, 0.2, 0.4, 1.0, 2.5 ]

# List of ionization parameters
logU_list = [-3, -2.5, -2]

# Name of template cloudy file
cloudy_input = osp.join('cloudy_slug', 'cloudy.in_grid_template')

# Name of the SLUG parameter file
slug_param = osp.join('param', 'cluster_cts_10myr.param')

# Loop over tracks and metallicities
for track, zrel in zip(track_list, zrel_list):

    # Edit the SLUG parameter file to match the desired track
    slug_param_tmp = '.'.join(slug_param.split('.')[:-1]) + \
                     '_' + '.'.join(track.split('.')[:-1]) + \
                     '.param'
    print "Writing SLUG parameter file "+slug_param_tmp
    fpin = open(slug_param, 'r')
    fpout = open(slug_param_tmp, 'w')
    for line in fpin:
        spl = line.split()
        if len(spl) > 0:
            if spl[0] == 'model_name':
                line = line[:-1] + '_' + '.'.join(track.split('.')[:-1])
                modname = line.split()[-1]
            elif spl[0] == 'tracks':
                line = 'tracks   lib/tracks/' + track
        fpout.write(line)
    fpin.close()
    fpout.close()

    # Edit the cloudy template to have the desired metallicity
    cloudy_input_tmp = cloudy_input + '_' + '.'.join(track.split('.')[:-1])
    print "Writing cloudy template file "+cloudy_input_tmp
    fpin = open(cloudy_input, 'r')
    fpout = open(cloudy_input_tmp, 'w')
    for line in fpin:
        fpout.write(line)
    fpout.write('metals and grains   {:f}'.format(zrel))
    fpin.close()
    fpout.close()

    # Run slug
    call([osp.join("bin", "slug"), slug_param_tmp])

    # Loop over ionization parameters
    for logU in logU_list:

        # Copy the slug files
        basename = modname+'_{:f}'.format(logU)
        shutil.copy(osp.join('output', modname+'_cluster_prop.fits'),
                    osp.join('output', basename + 
                             '_cluster_prop.fits'))
        shutil.copy(osp.join('output', modname+'_cluster_phot.fits'),
                    osp.join('output', basename + 
                             '_cluster_phot.fits'))
        shutil.copy(osp.join('output', modname+'_cluster_spec.fits'),
                    osp.join('output', basename + 
                             '_cluster_spec.fits'))

        # Use cloudy_slug to run cloudy
        call(["python", osp.join("cloudy_slug", "cloudy_slug.py"),
              basename, "--cloudytemplate", cloudy_input_tmp,
              "--clustermode", "--nodynamic",
              "--ionparam", "{:e}".format(10.**(logU)),
              "--nproc", "10",
              "--save", "--verbose"])

# Now read outputs
alldata = []

# Loop over tracks and metallicities
for track, zrel in zip(track_list, zrel_list):

    # Read the input file to get the model name
    slug_param_tmp = '.'.join(slug_param.split('.')[:-1]) + \
                     '_' + '.'.join(track.split('.')[:-1]) + \
                     '.param'
    fpin = open(slug_param, 'r')
    for line in fpin:
        spl = line.split()
        if len(spl) > 0:
            if spl[0] == 'model_name':
                line = line[:-1] + '_' + '.'.join(track.split('.')[:-1])
                modname = line.split()[-1]
                break
    fpin.close()

    # Loop over ionization parameters
    for logU in logU_list:

        # Read the cluster data
        basename = modname+'_{:f}'.format(logU)
        dataphot = read_cluster_phot(basename)
        datalines = read_cluster_cloudylines(basename)

        # Read the line array to get the exact wavelengths
        linearr_name = osp.join('cloudy_slug', 'cloudy_tmp_'+basename,
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
            if wl < 1e3:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + ' {:5.1f}A'.format(wl))
            elif wl < 1e4:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + '  {:4d}A'.format(int(wl)))
            elif wl < 1e5:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + ' {:5.3f}m'.format(wl/1e4))
            elif wl < 1e6:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + ' {:5.2f}m'.format(wl/1e4))
            else:
                linelabel.append(lab[0]+lab[1].lower()+lab[2:]
                                 + ' {:5.1f}m'.format(wl/1e4))
        linelabel = np.array(linelabel)

        # Identify entries that are not lines in the line array, and
        # throw these out
        idxkeep = []
        for i, lab in enumerate(linelabel):
            if lab in allwl:
                idxkeep.append(i)
        linelabel = linelabel[idxkeep]
        linelum = datalines.cloudy_linelum[:,idxkeep]

        # Identify H lines and throw these out
        idxkeep = []
        for i, lab in enumerate(linelabel):
            if lab[:2] != 'H ':
                idxkeep.append(i)
        linelabel = linelabel[idxkeep]
        linelum = datalines.cloudy_linelum[:,idxkeep]

        # Normalize all line luminosities to ionizing luminosity
        lumnorm = np.transpose(linelum)/dataphot.phot[:,0]

        # Identify lines that are always below 10^-20 / ionizing
        # photon in luminosity, and throw these out
        idxkeep = np.sum(lumnorm > 1e-20, axis=1) > 0
        linelabel = linelabel[idxkeep]
        lumnorm = lumnorm[idxkeep,:]

        # Get wavelengths for the remaining lines
        outwl=[]
        for lab in linelabel:
            outwl.append(allwl[lab])
        outwl = np.array(outwl)

        # Sort in order of decreasing luminosity at time 0
        idx = np.argsort(lumnorm[:,0])
        lumnorm = lumnorm[idx[::-1],:]
        outwl = outwl[idx[::-1]]
        linelabel = linelabel[idx[::-1]]

        # Append to master list
        outdict = { 
            'track' : track,
            'time' : dataphot.time,
            'logU' : logU,
            'wl' : outwl,
            'lum' : lumnorm,
            'label' : linelabel }
        alldata.append(outdict)

        #print "logU = {:f}, model = {:s}, Ha/QH0 = {:e}, nline = {:d}".format(
        #    logU, modname, 
        #    datalines.cloudy_linelum[0,53]/dataphot.phot[0,0],
        #    len(linelabel))

# Write out the metal lines file header
fp = open(osp.join(os.environ['SLUG_DIR'], 'lib', 'atomic', 'metallines.txt'), 'w')
fp.write('# Metal lines; this list was obtained by using SLUG to generate a\n')
fp.write('# fully-sampled spectrum for a Chabrier (2005) IMF at a range of ages,\n')
fp.write('# then using that spectrum in cloudy with constant density, background\n')
fp.write('# cosmic rays, expanding spehre geometry, and HII region abundances,\n')
fp.write('# scaled by  metallicity relative to Solar (which is taken to be 0.014\n')
fp.write('# for Geneva 2013 models, 0.020 for other models). The lines listed\n')
fp.write('# below are those for which the luminosity per ionizing photon\n')
fp.write('# injected is > 10^-20 erg/s/photon\n')
fp.write('#\n')
fp.write('# The first line below gives:\n')
fp.write('# NTRACK NTIME\n')
fp.write('# where NTRACK is the number of tracks and NTIME is the number of\n')
fp.write('# output times. The next line is NTIME entries of the form\n')
fp.write('# TIME_1 TIME_2 TIME_3 ...\n')
fp.write('# where TIME_N is the nth output time. The line after that is\n')
fp.write('# LOGU_1 LOGU_2 LOGU_3 ...\n')
fp.write('# where LOGU_N is the Nth value of the ionization parameter. This is\n')
fp.write('# followed by a series of NTRACKS * NU blocks, each of which begins:\n')
fp.write('# TRACK LOGU NLINES\n')
fp.write('# where TRACK is the track name, LOGU is the log of the ionization\n')
fp.write('# parameter, and NLINES is the number of lines for this block. This is\n')
fp.write('# followed by NLINES lines of the form\n')
fp.write('# SPEC     LAMBDA     (LUM/QH0)_1     (LUM/QH0)_2 ...\n')
fp.write('# where SPEC is the 11-character cloudy designation for the emitting\n')
fp.write('# species, LAMBDA is the line wavelength in Angstrom, and (LUM/QH0_N\n')
fp.write('# is the line luminosity per ionizing photon injected at the Nth time.\n')
fp.write('\n')

# Write out the data header
fp.write('{:d}  {:d}  {:d}\n'.format(len(track_list), len(alldata[0]['time']),
                                     len(logU_list)))
for t in alldata[0]['time']:
    fp.write('{:3.1e}   '.format(t))
fp.write('\n')
for logU in logU_list:
    fp.write('{:f}   '.format(logU))
fp.write('\n')

# Write out the data
for dat in alldata:
    fp.write('{:s}   {:f}   {:d}\n'.format(dat['track'], dat['logU'],
                                           len(dat['label'])))
    for label, wl, lum in zip(dat['label'], dat['wl'], dat['lum']):
        fp.write('{:s}   {:e}'.format(label, wl))
        for l in lum:
            fp.write('   {:e}'.format(l))
        fp.write('\n')

# Close the file
fp.close()
