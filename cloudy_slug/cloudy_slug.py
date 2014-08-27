"""
This is a script that takes spectra output by SLUG and calls cloudy on
them in order to calculate the resulting nebular emission.
"""

import argparse
import copy
import multiprocessing
import numpy as np
import os
import os.path as osp
from Queue import Queue
from scipy.constants import c
import subprocess
import sys
from threading import Thread
try:
    from slugpy import *    # If slugpy is already in our path
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")
alphaB = 2.59e-13

# Step 1: set up and read command line arguments
parser = argparse. \
         ArgumentParser(
             description="Script to run cloudy on slug outputs")
parser.add_argument("slug_model_name", 
                    help="name of the SLUG model output to be " +
                    "processed")
parser.add_argument("--slugpath", default=None, type=str,
                    help="path to the SLUG output data (default: "+
                    "checks cwd, then $SLUG_DIR/output)")
parser.add_argument("--cloudypath", default=None, type=str,
                    help="path to the cloudy executable (default: "+
                    "$CLOUDY_DIR/cloudy.exe)")
parser.add_argument("--cloudytemplate", default=None, type=str,
                    help="template cloudy input file (default: "+
                    "$SLUG_DIR/cloudy_slug/cloudy.in_template)")
parser.add_argument("--hden", nargs='+', metavar='nH', default=100,
                    help="hydrogen number densities (default: 100)")
parser.add_argument("-cm", "--clustermode", action='store_true',
                    default=False, help="run in cluster mode, where "+
                    "each cluster is a separate cloudy run "+
                    "(default: integrated mode, one cloudy run / "
                    "trial)")
parser.add_argument('-n', '--nproc', default=None, type=int,
                    help="number of cloudy processes (default: "+
                    "number of cores)")
parser.add_argument('-nl', '--nicelevel', default=0, type=int,
                    help="nice level of the cloudy processes " +
                    "(default: 0)")
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

# Step 3: read the SLUG output to be processed, and check that we have
# what we need
if (args.clustermode):
    data = read_cluster(args.slug_model_name, nofilterdata=True)
else:
    data = read_integrated(args.slug_model_name, nofilterdata=True)
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
freq = c/(data.wl*1e-10)          # Frequency in Hz
logfreq = np.log10(freq)          # Log frequency in Hz
basename = osp.basename(args.slug_model_name)


# Step 4: read the template cloudy input file template
fp = open(cloudytemplate, 'r')
tempfile = fp.read().split('\n')
fp.close()


# Step 5: queue up all the SLUG runs
slug_queue = Queue()
if (args.clustermode):
    for i in range(len(data.id)):
        slug_queue.put(i)
else:
    for i in range(data.spec.shape[-1]):
        for j in range(data.spec.shape[-2]):
            slug_queue.put((i, j))


# Step 6: define the worker function to do a cloudy run
def do_cloudy_run(thread_num, q):

    # Terminate when queue empties
    while not q.empty():

        # Fetch a task from the queue, get the data we need, and
        # generate a file name extension to go with it
        if args.clustermode:
            cluster_num = q.get()
            spec = data.spec[i,:]
            ext = "_n{:09d}".format(i)
            qH0 = data.phot[i,qH0idx]
            outstr = "launching cloudy on cluster {:d} of {:d}" \
                .format(i+1, len(data.id))
        else:
            i, j = q.get()
            spec = data.spec[:,j,i]
            ext = "_tr{:05d}_ti{:05d}".format(i,j)
            qH0 = data.phot[qH0idx,j,i]
            outstr = ("launching cloudy on trial {:d}/{:d}, " +
                      "time {:d}/{:d}").format(i+1, data.spec.shape[-1],
                                               j+1, data.spec.shape[-2])
            #print thread_num, i, j, ext

        # Write out the cloudy input file header, substituting a
        # custom name for OUTPUT_FILENAME
        cloudy_in_fname = osp.join(cwd, 
                                   'cloudy_tmp', 
                                   'cloudy.in'+'_{:05d}'.format(thread_num))
        fpout = open(cloudy_in_fname, 'w')
        radset = False
        hdenset = False
        outfiles = []
        for line in tempfile:
            linesplit = line.split()
            if len(linesplit) > 0:
                if (linesplit[0] == 'hden'):
                    hdenset = True
                    hden = 10.0**float(linesplit[1])
                elif (linesplit[0] == 'radius'):
                    radset = True
                if 'OUTPUT_FILENAME' in line:
                    fpout.write(
                        line.replace('OUTPUT_FILENAME',
                                     osp.join(cwd,
                                              'cloudy_tmp',
                                              basename+ext)) + '\n')
                else:
                    fpout.write(line+'\n')

        # Set H density and radius if setting them automatically
        if not hdenset:
            fpout.write("hden {:f}\n".format(np.log10(args.hden)))
            hden = args.hden
        if not radset:
            rstrom = (3.0*qH0/(4.0*np.pi*alphaB*hden**2))**(1./3.)
            fpout.write("radius {:f}\n".format(np.log10(rstrom)-3.0))

        # Write the ionizing or bolometric luminosity to the cloudy
        # input file
        if qH0idx != -1:
            fpout.write("Q(H) = {:f}\n".format(np.log10(qH0)))
        else:
            fpout.write("luminosity solar {:f}\n".
                        format(np.log10(lbol)))

        # Write the spectral shape into the cloudy input file,
        # prepending and appending low values outside the range
        # covered by SLUG to make cloudy happy
        logL_nu = np.log10(spec*c/freq**2)
        fpout.write("interpolate")
        fpout.write(" ({:f} {:f})".format(7.51, 
                                          np.amin(logL_nu)-4))
        fpout.write(" ({:f} {:f})".format(logfreq[-1]-0.01, 
                                          np.amin(logL_nu)-4))
        for i in range(len(logL_nu)):
            if i % 4 == 0:
                fpout.write("\ncontinue")
            fpout.write(" ({:f} {:f})".format(logfreq[-i-1],
                                              logL_nu[-i-1]))
        fpout.write("\ncontinue ({:f} {:f})".
                    format(logfreq[0]+0.01, 
                           np.amin(logL_nu)-4))
        fpout.write(" ({:f} {:f})\n".format(22.4, 
                                            np.amin(logL_nu)-4))

        # Close cloudy input file
        fpout.close()

        # If verbose, print status
        if args.verbose:
            print("thread {:d}: ".format(thread_num) + outstr)

        # Launch the cloudy process
        cmd = cloudypath + " < " + cloudy_in_fname + \
              " > " + osp.join(cwd, 'cloudy_tmp', 'cloudy.out'+ext)
        if args.nicelevel > 0:
            cmd = "nice -n " + str(args.nicelevel) + " " + cmd
        #print cmd
        proc = subprocess.Popen(cmd, shell=True)

        # Block until the cloudy process completes, then declare that
        # we're done
        proc.wait()
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


# Step 8: 
