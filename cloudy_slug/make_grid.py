# This script generates the grid that SLUG uses for its "quick and
# dirty" estimate of metal line luminosities and HII region
# temperatures.

import os
import os.path as osp
import sys
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
zrel_list = [
    # MIST non-rotating
    10.**-4.0,
    10.**-3.5, 10.**-3.0, 10.**-2.5, 10.**-2.0, 10.**-1.75, 10.**-1.5, 10.**-1.25, 10.**-1.0, 10.**-0.75, 10.**-0.5,
    10.**-0.25, 10.**0.0, 10.**0.25, 10.**0.5,
    # MIST rotating
    10.**-4.0, 10.**-3.5, 10.**-3.0, 10.**-2.5, 10.**-2.0, 10.**-1.75, 10.**-1.5, 10.**-1.25, 10.**-1.0, 10.**-0.75, 10.**-0.5,
    10.**-0.25, 10.**0.0, 10.**0.25, 10.**0.5,
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
cloudy_input = osp.join(os.environ['SLUG_DIR'],
                        'cloudy_slug', 'cloudy.in_grid_template')

# Name of the SLUG parameter files
slug_param = osp.join(os.environ['SLUG_DIR'], 'param',
                      'cluster_cts_10myr.param')
slug_param_gal = osp.join(os.environ['SLUG_DIR'],
                          'param', 'galaxy_cts_10myr.param')

# Loop over tracks and metallicities
for track, zrel in zip(track_list, zrel_list):

    # Edit the SLUG parameter files to match the desired track
    slug_param_tmp = '.'.join(slug_param.split('.')[:-1]) + \
                     '_' + '.'.join(
                         osp.basename(track).split('.')[:-1]) + \
                     '.param'
    print("Writing SLUG parameter file "+slug_param_tmp)
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
    print("Writing SLUG parameter file "+slug_param_tmp_gal)
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

    # Edit the cloudy template to have the desired metallicity
    cloudy_input_tmp = cloudy_input + '_' + '.'.join(
        osp.basename(track).split('.')[:-1])
    print("Writing cloudy template file "+cloudy_input_tmp)
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
    if not osp.isfile(osp.join("nebular_grid",
                               modname_gal+"_integrated_spec.fits")):
        call([osp.join(os.environ['SLUG_DIR'],
                       "bin", "slug"), slug_param_tmp_gal])
    else:
        print("Found existing model "+modname_gal+
              ", skipping slug run")

    # Run slug on the cluster case
    if not osp.isfile(osp.join("nebular_grid",
                               modname+"_cluster_spec.fits")):
        call([osp.join(os.environ['SLUG_DIR'],
                       "bin", "slug"), slug_param_tmp])
    else:
        print("Found existing model "+modname+
              ", skipping slug run")

    # Loop over ionization parameters
    for logU in logU_list:

        # Make names
        basename = modname+'_{:f}'.format(logU)
        basename_gal = modname_gal+'_{:f}'.format(logU)

        # Copy the slug files
        shutil.copy(osp.join('nebular_grid', modname+'_cluster_prop.fits'),
                    osp.join('nebular_grid', basename + 
                             '_cluster_prop.fits'))
        shutil.copy(osp.join('nebular_grid', modname+'_cluster_phot.fits'),
                    osp.join('nebular_grid', basename + 
                             '_cluster_phot.fits'))
        shutil.copy(osp.join('nebular_grid', modname+'_cluster_spec.fits'),
                    osp.join('nebular_grid', basename + 
                             '_cluster_spec.fits'))
        shutil.copy(osp.join('nebular_grid', modname_gal+'_integrated_prop.fits'),
                    osp.join('nebular_grid', basename_gal + 
                             '_integrated_prop.fits'))
        shutil.copy(osp.join('nebular_grid', modname_gal+'_integrated_phot.fits'),
                    osp.join('nebular_grid', basename_gal + 
                             '_integrated_phot.fits'))
        shutil.copy(osp.join('nebular_grid', modname_gal+'_integrated_spec.fits'),
                    osp.join('nebular_grid', basename_gal + 
                             '_integrated_spec.fits'))

        # Use cloudy_slug to run cloudy on the integrated case; skip if already done
        if not osp.isfile(
                osp.join("nebular_grid",
                         basename_gal+"_integrated_cloudylines.fits")):
            cmdstr = sys.executable + " " + \
                     osp.join(os.environ['SLUG_DIR'], 
                              "cloudy_slug", "cloudy_slug.py") + " " + \
                     osp.join("nebular_grid", basename_gal) + " " + \
                     "--cloudytemplate " + cloudy_input_tmp + " " + \
                     "--ionparam {:e}".format(10.**(logU)) + " " + \
                     "--hden {:e}".format(100.0) + " " + \
                     "--save --warnfail --verbose --fix_quantity nII"
            print(cmdstr)
            call(cmdstr, shell=True)
        else:
            print("Found cloudy files for "+
                  basename_gal+"; skipping...")

        # Use cloudy_slug to run cloudy on the cluster case; skip if already done
        if not osp.isfile(
                osp.join("nebular_grid",
                         basename+"_cluster_cloudylines.fits")):
            cmdstr = sys.executable + " " + \
                     osp.join(os.environ['SLUG_DIR'],
                              "cloudy_slug", "cloudy_slug.py") + " " + \
                     osp.join("nebular_grid", basename) + " " + \
                     "--cloudytemplate " + cloudy_input_tmp + " " + \
                     "--ionparam {:e}".format(10.**(logU)) + " " + \
                     "--hden {:e}".format(100.0) + " " + \
                     "--clustermode --fix_quantity nII " + \
                     "--save --warnfail --verbose"
            print(cmdstr)
            call(cmdstr, shell=True)
        else:
            print("Found cloudy files for "+
                  basename+"; skipping...")


