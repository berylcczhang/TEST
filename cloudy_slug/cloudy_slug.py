"""
This is a script that takes spectra output by SLUG and calls cloudy on
them in order to calculate the resulting nebular emission.
"""

import argparse
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
from collections import namedtuple
import copy
import multiprocessing
import numpy as np
import os
import os.path as osp
try:
    from Queue import Queue    # python 2.x
except ImportError:
    from queue import Queue    # python 3.x
from scipy.constants import c
from scipy.constants import k as kB
import subprocess
import sys
from threading import Thread
from time import sleep
try:
    from slugpy import *    # If slugpy is already in our path
    from slugpy.cloudy import *
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        from slugpy.cloudy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Step 1: set up and read command line arguments
parser = argparse. \
         ArgumentParser(
             description="Script to run cloudy on slug outputs")

# Positional arguments
parser.add_argument("slug_model_name", 
                    help="name of the SLUG model output to be " +
                    "processed")
parser.add_argument("start_spec", nargs="?", type=int, default=0,
                    help="starting cluster or trial number (default: 0)")
parser.add_argument("end_spec", nargs="?", type=int, default=-1,
                    help="ending cluster or trial number " + 
                    "(default: last spectrum)")

# Optional arguments
parser.add_argument('-a', '--agemax', default=4, type=float,
                    help="maximum cluster age in Myr for which to " +
                    "compute nebular emission; only used in clustermode " +
                    "(default: 4 Myr)")
parser.add_argument("--cloudypath", default=None, type=str,
                    help="path to the cloudy executable (default: "+
                    "$CLOUDY_DIR/cloudy.exe)")
parser.add_argument("--cloudytemplate", default=None, type=str,
                    help="template cloudy input file (default: "+
                    "$SLUG_DIR/cloudy_slug/cloudy.in_template)")
parser.add_argument("-cm", "--clustermode", action='store_true',
                    default=False, help="run in cluster mode, where "+
                    "each cluster is a separate cloudy run "+
                    "(default: integrated mode, one cloudy run / "
                    "trial)")
parser.add_argument('-nl', '--nicelevel', default=0, type=int,
                    help="nice level of the cloudy processes " +
                    "(default: 0)")
parser.add_argument('-n', '--nproc', default=None, type=int,
                    help="number of cloudy processes (default: "+
                    "number of cores)")
parser.add_argument('-s', '--save', default=False,
                    action='store_true', help='save full cloudy ' +
                    'output (default: delete after extracting data)')
parser.add_argument("--slugpath", default=None, type=str,
                    help="path to the SLUG output data (default: "+
                    "checks cwd, then $SLUG_DIR/output)")
parser.add_argument('-v', '--verbose', action='store_true',
                    default=False, help="produce verbose output")
args = parser.parse_args()
cwd = osp.dirname(osp.realpath(__file__))


# Step 2: set some defaults
if args.cloudypath is None:
    if 'CLOUDY_DIR' in os.environ:
        cloudypath = osp.join(os.environ['CLOUDY_DIR'], 'cloudy.exe')
    else:
        cloudypath = 'cloudy.exe'
else:
    cloudypath = args.cloudypath
if args.nproc is None:
    nproc = multiprocessing.cpu_count()
else:
    nproc = args.nproc
if args.cloudytemplate is None:
    if 'SLUG_DIR' in os.environ:
        cloudytemplate = osp.join(os.environ['SLUG_DIR'],
                                  'cloudy_slug',
                                  'cloudy.in_template')
    else:
        cloudytemplate = osp.join('cloudy_slug',
                                  'cloudy.in_template')
else:
    cloudytemplate = args.cloudytemplate

if args.slugpath is None:
    #set input directory 
    if 'SLUG_DIR' in os.environ:
        slugdir=osp.join(os.environ['SLUG_DIR'],'output/')
    else:
        slugdir='output/'
    #use cwd as output
    outpath='./'
else:
    slugdir=args.slugpath
    #use selected output dir for  output
    outpath=slugdir

# Step 3: read the SLUG output to be processed, and check that we have
# what we need
file_info = {}
if (args.clustermode):
    data = read_cluster(args.slug_model_name, read_info=file_info, output_dir=slugdir)
else:
    data = read_integrated(args.slug_model_name, read_info=file_info, output_dir=slugdir)
valid = True
if 'spec' not in data._fields:
    valid = False
if 'filter_names' not in data._fields:
    valid = False
else:
    if 'QH0' not in data.filter_names:
        valid = False
    elif 'QH0' in data.filter_names:
        qH0idx = data.filter_names.index('QH0')
if not valid:
    raise IOError("cloudy_slug: error: input slug data must " +
                  "contain spectra and ionizing luminosity")
if args.clustermode and 'form_time' not in data._fields:
    raise IOError("cloudy_slug: error: input slug data must " +
                  "contain cluster physical properties")
freq = c/(data.wl*1e-10)          # Frequency in Hz
logfreq = np.log10(freq)          # Log frequency in Hz
basename = osp.basename(args.slug_model_name)
if args.clustermode:
    if args.end_spec != -1:
        end_spec = min(args.end_spec, len(data.id))
    else:
        end_spec = len(data.id)
else:
    if args.end_spec != -1:
        end_spec = min(args.end_spec, data.spec.shape[-1])
    else:
        end_spec = data.spec.shape[-1]

# Figure out the photometric system we're using
if 'erg/s/Hz' in data.filter_units:
    photsystem = 'L_nu'
elif 'erg/s/A' in data.filter_units:
    photsystem = 'L_lambda'
elif 'AB mag' in data.filter_units:
    photsystem = 'AB'
elif 'ST mag' in data.filter_units:
    photsystem = 'STMAG'
elif 'Vega mag' in data.filter_units:
    photsystem = 'Vega'
else:
    photsystem = 'L_nu'

# Number of filters
nfilt = len(data.filter_names)


# Step 4: read the template cloudy input file template, and set up
# storage for what we'll be computing
fp = open(cloudytemplate, 'r')
tempfile = fp.read().split('\n')
fp.close()
compute_continuum = False
compute_lines = False
for line in tempfile:
    if 'save' in line and 'continuum' in line:
        compute_continuum = True
    elif 'save' in line and 'lines' and 'list' in line:
        compute_lines = True
if compute_continuum:
    cloudywl = []
    cloudyspec = []
    cloudyphot = []
    # Create dummy holders for cloudy spectra
    if args.clustermode:
        for i in range(end_spec-args.start_spec):
            cloudywl.append(None)
            cloudyspec.append(None)
            cloudyphot.append(None)
    else:
        for j in range(data.spec.shape[-2]):
            cloudywl.append([])
            cloudyspec.append([])
            cloudyphot.append([])
            for i in range(args.start_spec, end_spec):
                cloudywl[j].append(None)
                cloudyspec[j].append(None)
                cloudyphot[j].append(None)
if compute_lines:
    linelist = None
    linewl = None
    linelum = []
    # Create dummy holders for cloudy lines
    if args.clustermode:
        for i in range(end_spec-args.start_spec):
            linelum.append(None)
    else:
        for j in range(data.spec.shape[-2]):
            linelum.append([])
            for i in range(args.start_spec, end_spec):
                linelum[j].append(None)


# Step 5: queue up the SLUG runs
slug_queue = Queue()
if args.clustermode:
    for i in range(args.start_spec, end_spec):
        slug_queue.put(i)
else:
    for i in range(args.start_spec, end_spec):
        for j in range(data.spec.shape[-2]):
            slug_queue.put((i, j))

# Step 6: define the worker function to do a cloudy run
def do_cloudy_run(thread_num, q):

    # Declare global variables
    global compute_continuum
    global compute_lines
    global data
    if compute_continuum:
        global cloudywl
        global cloudyspec
        global cloudyphot
    if compute_lines:
        global linelist
        global linewl
        global linelum

    # Terminate when queue empties
    while not q.empty():

        # Fetch a task from the queue, get the data we need, and
        # generate a file name extension to go with it
        if args.clustermode:

            # Cluster mode
            cluster_num = q.get()

            # Check if this cluster is below our age maximum; if not,
            # skip it
            if data.time[cluster_num]-data.form_time[cluster_num] > \
               args.agemax*1e6*365.25*24.*3600.:
                q.task_done()
                continue
            spec = data.spec[cluster_num,:]
            ext = "_n{0:09d}".format(cluster_num)
            qH0 = data.phot[cluster_num,qH0idx]
            outstr = "launching cloudy on cluster {0:d} of {1:d}" \
                .format(cluster_num+1, len(data.id))

        else:

            # Integrated mode
            trial, time = q.get()
            spec = data.spec[:, time, trial]
            ext = "_tr{0:05d}_ti{1:05d}".format(trial,time)
            qH0 = data.phot[qH0idx,time,trial]
            outstr = ("launching cloudy on trial {0:d}/{1:d}, " +
                      "time {2:d}/{3:d}").format(trial+1, data.spec.shape[-1],
                                               time+1, data.spec.shape[-2])

        # Write out the cloudy input file header, substituting a
        # custom name for OUTPUT_FILENAME
        cloudy_in_fname = osp.join(cwd, 
                                   'cloudy_tmp', 
                                   'cloudy.in'+'_{0:05d}'.format(thread_num))
        fpout = open(cloudy_in_fname, 'w')
        radset = False
        for line in tempfile:
            linesplit = line.split()
            if len(linesplit) > 0:
                if (linesplit[0] == 'hden'):
                    hden = 10.0**float(linesplit[1])
                    if not args.clustermode:
                        fpout.write(line+'\n')
                elif (linesplit[0] == 'radius'):
                    radset = True
                    fpout.write(line+'\n')
                elif 'OUTPUT_FILENAME' in line:
                    newline \
                        = line.replace('OUTPUT_FILENAME',
                                       osp.join(cwd, 'cloudy_tmp',
                                                basename+ext))
                    fpout.write(newline + '\n')
                    if 'continuum' in newline:
                        lquote = newline.find('"')
                        rquote = newline.rfind('"')
                        continuum_file = newline[lquote+1:rquote]
                    elif 'line list' in newline or 'linelist' in newline:
                        lquote = newline.find('"')
                        rquote = newline[lquote+1:].find('"')
                        lines_file = newline[lquote+1:lquote+1+rquote]
                else:
                    fpout.write(line+'\n')

        # In cluster mode, compute internal radius, starting density
        # self-consistently; in integrated mode, set inner radius
        # self-consistently
        if not args.clustermode:
            if not radset:
                alphaB = 2.59e-13    # Case B recombination coefficient
                rstrom = (3.0*qH0/(4.4*np.pi*alphaB*hden**2))**(1./3.)
                r0 = rstrom/1e3
                fpout.write("radius {0:f}\n".format(np.log10(r0)))
        else:
            # Get current age of this cluster
            age = data.time[0] - data.form_time[0]
            # Get characteristic radius and time; all quantities
            # defined as in Krumholz & Matnzer (2009)
            alphaB = 2.59e-13    # Case B recombination coefficient
            mu = 2.34e-24        # Mean mass per H nucleus, in g
            eps0 = 2.179e-11     # H ionization potential, in erg
            TII = 1e4            # HII region temp, in K
            psi = 3.2            # Mean photon energy / eps0
            ft = 2.0             # Trapping factor
            phi = 0.73           # Dust absorption fraction
            rch = alphaB/(12.0*np.pi*phi) * \
                  (eps0/(2.2*kB*1e7*TII))**2 * ft**2 * \
                  psi**2 * qH0 / (c*1e2)**2
            tch = (4.0*np.pi * mu*hden*c * rch**4 /
                   (3.0*ft*qH0*psi*eps0))**0.5
            # Get xIIgas, xIIrad
            tau = age*365.25*24.*3600./tch
            xIIrad = (2.0*tau**2)**0.25
            xIIgas = (49.0*tau**2/36.0)**(2.0/7.0)
            # Get outer radius, inner radius, density
            r = rch*(xIIrad**3.5 + xIIgas**3.5)**(2.0/7.0)
            r0 = r/1e3
            nH = (3.0*qH0 / (4.0*np.pi*alphaB*r**3))**0.5
            fpout.write("hden {0:f}\n".format(np.log10(nH)))
            if not radset:
                fpout.write("radius {0:f}\n".format(np.log10(r0)))

        # Write the ionizing luminosity to the cloudy input file
        fpout.write("Q(H) = {0:f}\n".format(np.log10(qH0)))

        # Write the spectral shape into the cloudy input file,
        # prepending and appending low values outside the range
        # covered by SLUG to make cloudy happy. Make sure to handle
        # zeros gracefully.
        specclean = np.copy(spec)
        specclean[spec == 0.0] = np.amin(spec[spec > 0])*1e-4
        logL_nu = np.log10(specclean*c/freq**2)
        fpout.write("interpolate")
        fpout.write(" ({0:f} {1:f})".format(7.51, 
                                          np.amin(logL_nu)-4))
        fpout.write(" ({0:f} {1:f})".format(logfreq[-1]-0.01, 
                                          np.amin(logL_nu)-4))
        for i in range(len(logL_nu)):
            if i % 4 == 0:
                fpout.write("\ncontinue")
            fpout.write(" ({0:f} {1:f})".format(logfreq[-i-1],
                                              logL_nu[-i-1]))
        fpout.write("\ncontinue ({0:f} {1:f})".
                    format(logfreq[0]+0.01, 
                           np.amin(logL_nu)-4))
        fpout.write(" ({0:f} {1:f})\n".format(22.4, 
                                            np.amin(logL_nu)-4))

        # Close cloudy input file
        fpout.close()

        # If verbose, print status
        if args.verbose:
            print("thread {0:d}: ".format(thread_num+1) + outstr)

        # Launch the cloudy process and wait for it to complete
        cloudy_out_fname = osp.join(cwd, 'cloudy_tmp', 'cloudy.out'+ext)
        cmd = cloudypath + " < " + cloudy_in_fname + \
              " > " + cloudy_out_fname
        if args.nicelevel > 0:
            cmd = "nice -n " + str(args.nicelevel) + " " + cmd
        proc = subprocess.call(cmd, shell=True)

        # Read and store the cloudy continuum output
        if continuum_file is not None:
            while os.stat(continuum_file).st_size == 0:
                sleep(2)
            cdata = read_cloudy_continuum(continuum_file, r0=r0)
            if args.clustermode:
                cloudywl[cluster_num] = cdata.wl
                cloudyspec[cluster_num] = cdata.L_lambda
            else:
                cloudywl[time][trial] = cdata.wl
                cloudyspec[time][trial] = cdata.L_lambda

        # Read and store the cloudy line luminosity output
        if lines_file is not None:
            while os.stat(lines_file).st_size == 0:
                sleep(2)
            ldata = read_cloudy_linelist(lines_file)
            if linelist is None:
                linelist = ldata.label
            if linewl is None:
                linewl = ldata.wl
            if args.clustermode:
                linelum[cluster_num] = ldata.lum
            else:
                linelum[time][trial] = ldata.lum

        # Compute photometry from the cloudy data
        if continuum_file is not None:
            trans_phot \
                = compute_photometry(cdata.wl, cdata.L_lambda[1,:],
                                     data.filter_names, 
                                     photsystem=photsystem,
                                     filter_wl=data.filter_wl,
                                     filter_response=data.filter_response,
                                     filter_beta=data.filter_beta,
                                     filter_wl_c=data.filter_wl_c)
            emit_phot \
                = compute_photometry(cdata.wl, cdata.L_lambda[2,:],
                                     data.filter_names, 
                                     photsystem=photsystem,
                                     filter_wl=data.filter_wl,
                                     filter_response=data.filter_response,
                                     filter_beta=data.filter_beta,
                                     filter_wl_c=data.filter_wl_c)
            trans_emit_phot \
                = compute_photometry(cdata.wl, cdata.L_lambda[3,:],
                                     data.filter_names, 
                                     photsystem=photsystem,
                                     filter_wl=data.filter_wl,
                                     filter_response=data.filter_response,
                                     filter_beta=data.filter_beta,
                                     filter_wl_c=data.filter_wl_c)
            phot = np.array([trans_phot, emit_phot,
                             trans_emit_phot])
            if args.clustermode:
                cloudyphot[cluster_num] = phot
            else:
                cloudyphot[time][trial] = phot

        # Clean up the cloudy output unless requested to keep it
        if not args.save:
            if continuum_file is not None:
                os.remove(continuum_file)
            if lines_file is not None:
                os.remove(lines_file)
            os.remove(cloudy_in_fname)
            os.remove(cloudy_out_fname)

        # Declare that we're done
        q.task_done()


# Step 7: start a set of threads to do the job
try: 
    os.mkdir(osp.join(cwd, 'cloudy_tmp'))  # Temporary working directory
except OSError: pass           # Probably failed because dir exists
if __name__ == '__main__' :
    for i in range(nproc):
        p = Thread(target=do_cloudy_run, args=(i, slug_queue))
        p.start()
slug_queue.join()


# Step 8: convert the continuum spectra to a namedtuple in the same
# format as slug spectra. Note that we have to pad the continuum
# spectra we get from cloudy to make them all the same size, so
# that the can be converted to arrays; the continuum spectra are
# assumed to follow standard cloudy output format, so that they all
# have the same wavelength spacing and the same maximum wavelength,
# but different minimum wavelengths. Thus we are padding the
# beginnings of the arrays.
if compute_continuum:

    # Get the maximum length wavelength array
    cloudywl_max = np.zeros(0)
    if args.clustermode:
        for i in range(len(cloudywl)):
            if cloudy_wl[i] is not None:
                if cloudywl[i].shape[0] > cloudywl_max.shape[0]:
                    cloudywl_max = cloudywl[i]
    else:
        for i in range(len(cloudywl)):
            for j in range(len(cloudywl[0])):
                if cloudywl[i][j].shape[0] > cloudywl_max.shape[0]:
                    cloudywl_max = cloudywl[i][j]

    # Now loop over stored spectra, padding array beginnings
    if args.clustermode:
        for i in range(len(cloudywl)):
            if cloudy_wl[i] is not None:
                offset = len(cloudywl_max) - len(cloudywl[i])
                if offset > 0:
                    cloudyspec[i] \
                        = np.insert(cloudyspec[i], 0,
                                    np.zeros((offset,4)), axis=1)
            else:
                cloudyspec[i] = np.zeros(len(cloudywl_max))
    else:
        for i in range(len(cloudywl)):
            for j in range(len(cloudywl[0])):
                offset = len(cloudywl_max) - len(cloudywl[i][j])
                if offset > 0:
                    cloudyspec[i][j] \
                        = np.insert(cloudyspec[i][j], 0,
                                    np.zeros((offset,4)), axis=1)

    # Now that we've made the spectra the same size, turn them into an
    # array
    cloudyspec = np.array(cloudyspec)

    # Final step: make a namedtuple to hold the data
    if args.clustermode:
        cloudyspec_type = namedtuple('cluster_cloudyspec',
                                     ['id', 'trial', 'time', 
                                      'cloudy_wl', 'cloudy_inc', 
                                      'cloudy_trans', 'cloudy_emit', 
                                      'cloudy_trans_emit'])
        cloudyspec_data \
            = cloudyspec_type(data.id, data.trial, data.time, cloudywl_max,
                              cloudyspec[:,0,:], cloudyspec[:,1,:],
                              cloudyspec[:,2,:], cloudyspec[:,3,:])
    else:
        cloudyspec_type = namedtuple('integrated_cloudyspec',
                                     ['time', 'cloudy_wl', 'cloudy_inc', 
                                      'cloudy_trans', 'cloudy_emit', 
                                      'cloudy_trans_emit'])
        cloudyspec_data \
            = cloudyspec_type(data.time, cloudywl_max,
                              np.transpose(cloudyspec[:,:,0,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,1,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,2,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,3,:], (2,0,1)))

# Step 9: write the cloudy spectra to file
if compute_continuum:
    if args.clustermode:
        write_cluster_cloudyspec(cloudyspec_data, osp.join(outpath,args.slug_model_name),
                                 file_info['format'])
    else:
        write_integrated_cloudyspec(cloudyspec_data,
                                    osp.join(outpath,args.slug_model_name),
                                    file_info['format'])

# Step 10: write the line data to file
if compute_lines:
    if args.clustermode:
 
        # Set clusters we skipped to have line luminosities of zero
        nline = 1
        for i in range(len(linelum)):
            if linelum[i] is not None:
                nline = linelum[i].shape[0]
                break
        for i in range(len(linelum)):
            if linelum[i] is None:
                linelum[i] = np.zeros(nline)

        # Write line data
        linelum = np.array(linelum)
        cloudylines_type = namedtuple('cluster_cloudylines',
                                      ['id', 'trial', 'time', 
                                       'cloudy_linelist',
                                       'cloudy_linewl',
                                       'cloudy_linelum'])
        cloudylines = cloudylines_type(data.id, data.trial, data.time, 
                                       linelist, linewl, linelum)
        write_cluster_cloudylines(cloudylines, osp.join(outpath,args.slug_model_name),
                                  file_info['format'])
    else:
        linelum = np.array(linelum)
        cloudylines_type = namedtuple('cluster_cloudylines',
                                      ['time', 'cloudy_linelabel',
                                       'cloudy_linewl',
                                       'cloudy_linelum'])
        cloudylines = cloudylines_type(data.time, linelist,
                                       linewl, 
                                       np.transpose(linelum, (2,0,1)))
        write_integrated_cloudylines(cloudylines, osp.join(outpath,args.slug_model_name),
                                     file_info['format'])

# Step 11: write photometry to file
if compute_continuum:

    # Cluster or integrated mode
    if args.clustermode:

        # Set clusters we skipped to have photometric values of either
        # 0 (for non-magnitude systems) or +infinity (for magnitude
        # systems)
        nfilter = len(data.filter_names)
        for i in range(len(cloudyphot)):
            if cloudyphot[i] is None:
                cloudyphot[i] = np.zeros((3,nfilter))
                for j in range(nfilter):
                    if 'mag' in data.filter_units[j]:
                        cloudyphot[i][:,j] = np.inf

        # Convert to array
        cloudyphot = np.array(cloudyphot)

        # Build namedtuple
        cloudyphot_type = namedtuple('cluster_cloudyphot',
                                     ['id', 'trial', 'time', 
                                      'cloudy_filter_names', 
                                      'cloudy_filter_units',
                                      'cloudy_filter_wl_eff', 
                                      'cloudy_filter_wl',
                                      'cloudy_filter_response',
                                      'cloudy_filter_beta',
                                      'cloudy_filter_wl_c',
                                      'cloudy_phot_trans',
                                      'cloudy_phot_emit', 
                                      'cloudy_phot_trans_emit'])
        cloudyphot_data \
            = cloudyphot_type(data.id, data.trial, data.time,
                              data.filter_names, data.filter_units,
                              data.filter_wl_eff, data.filter_wl,
                              data.filter_response, data.filter_beta,
                              data.filter_wl_c, cloudyphot[:,0,:],
                              cloudyphot[:,1,:], cloudyphot[:,2,:])

        # Write
        write_cluster_cloudyphot(cloudyphot_data, osp.join(outpath,args.slug_model_name),
                                 file_info['format'])

    else:

        # Convert to array
        cloudyphot = np.array(cloudyphot)

        # Build namedtuple
        cloudyphot_type = namedtuple('integrated_cloudyphot',
                                     ['time', 
                                      'cloudy_filter_names', 
                                      'cloudy_filter_units',
                                      'cloudy_filter_wl_eff', 
                                      'cloudy_filter_wl',
                                      'cloudy_filter_response',
                                      'cloudy_filter_beta',
                                      'cloudy_filter_wl_c',
                                      'cloudy_phot_trans',
                                      'cloudy_phot_emit', 
                                      'cloudy_phot_trans_emit'])
        cloudyphot_data \
            = cloudyphot_type(data.time, 
                              data.filter_names, data.filter_units,
                              data.filter_wl_eff, data.filter_wl,
                              data.filter_response, data.filter_beta,
                              data.filter_wl_c,
                              np.transpose(cloudyphot[:,:,0,:], (2,0,1)),
                              np.transpose(cloudyphot[:,:,1,:], (2,0,1)),
                              np.transpose(cloudyphot[:,:,2,:], (2,0,1)))

        # Write
        write_integrated_cloudyphot(cloudyphot_data, osp.join(outpath,args.slug_model_name),
                                    file_info['format'])


# Step 12: final cleanup
if not args.save:
    try:
        os.rmdir(osp.join(cwd, 'cloudy_tmp'))
    except OSError:
        warnings.warn("unable to clean up temporary directory "+
                      osp.join(cwd, 'cloudy_tmp'))
