# This script generates the grid that SLUG uses for its "quick and
# dirty" estimate of metal line luminosities and HII region
# temperatures.

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
    # Geneva 2013 non-rotating models
    'Z0140v00.txt', 'Z0020v00.txt', 
    # Padova with TP-AGB models
    'modp0004.dat', 'modp004.dat', 'modp008.dat', 'modp020.dat',
    'modp050.dat',
    # Geneva 2014 rotating models
    'Z0020v40.txt', 'Z0140v40.txt',
    # Geneva pre-2013, standard mass loss
    'modc001.dat', 'modc004.dat', 'modc008.dat', 'modc020.dat', 
    'modc040.dat',
    # Geneva pre-2013, enhanced mass los
    'mode001.dat', 'mode004.dat', 'mode008.dat', 'mode020.dat', 
    'mode040.dat',
    # Padova w/out TP-AGB
    'mods0004.dat', 'mods004.dat', 'mods008.dat', 'mods020.dat',
    'mods050.dat'
]
zrel_list = [
    # Geneva 2013 non-rotating models
    1.0, 1.0/7.0, 
    # Padova with TP-AGB models
    0.02, 0.2, 0.4, 1.0, 2.5,
    # Geneva 2014 rotating models
    1.0, 1.0/7.0, 
    # Geneva pre-2013, standard mass loss
    0.05, 0.2, 0.4, 1.0, 2.0,
    # Geneva pre-2013, enhanced mass los
    0.05, 0.2, 0.4, 1.0, 2.0,
    # Padova w/out TP-AGB
    0.02, 0.2, 0.4, 1.0, 2.5
 ]

# List of ionization parameters
logU_list = [-3, -2.5, -2]

# Name of template cloudy file
cloudy_input = osp.join('cloudy_slug', 'cloudy.in_grid_template')

# Name of the SLUG parameter files
slug_param = osp.join('param', 'cluster_cts_10myr.param')
slug_param_gal = osp.join('param', 'galaxy_cts_10myr.param')

# Loop over tracks and metallicities
#for i in range(0):
for track, zrel in zip(track_list, zrel_list):

    # Edit the SLUG parameter files to match the desired track
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
    slug_param_tmp_gal = '.'.join(slug_param_gal.split('.')[:-1]) + \
                     '_' + '.'.join(track.split('.')[:-1]) + \
                     '.param'
    print "Writing SLUG parameter file "+slug_param_tmp_gal
    fpin = open(slug_param_gal, 'r')
    fpout = open(slug_param_tmp_gal, 'w')
    for line in fpin:
        spl = line.split()
        if len(spl) > 0:
            if spl[0] == 'model_name':
                line = line[:-1] + '_' + '.'.join(track.split('.')[:-1])
                modname_gal = line.split()[-1]
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
        if line[-1] != '\n':
            fpout.write('\n')
    fpout.write('metals and grains   {:f}'.format(zrel))
    fpin.close()
    fpout.close()

    # Run slug on the galaxy case
    call([osp.join("bin", "slug"), slug_param_tmp_gal])

    # Run slug on the cluster case
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
        basename_gal = modname_gal+'_{:f}'.format(logU)
        shutil.copy(osp.join('output', modname_gal+'_integrated_prop.fits'),
                    osp.join('output', basename_gal + 
                             '_integrated_prop.fits'))
        shutil.copy(osp.join('output', modname_gal+'_integrated_phot.fits'),
                    osp.join('output', basename_gal + 
                             '_integrated_phot.fits'))
        shutil.copy(osp.join('output', modname_gal+'_integrated_spec.fits'),
                    osp.join('output', basename_gal + 
                             '_integrated_spec.fits'))

        # Use cloudy_slug to run cloudy on the integrated case
        call(["python", osp.join("cloudy_slug", "cloudy_slug.py"),
              basename_gal, "--cloudytemplate", cloudy_input_tmp,
              "--ionparam", "{:e}".format(10.**(logU)),
              "--save", "--verbose"])

        # Use cloudy_slug to run cloudy on the cluster case
        call(["python", osp.join("cloudy_slug", "cloudy_slug.py"),
              basename, "--cloudytemplate", cloudy_input_tmp,
              "--clustermode", "--nodynamic",
              "--ionparam", "{:e}".format(10.**(logU)),
              "--nproc", "10",
              "--save", "--verbose"])

