#!/usr/bin/env python

"""
This is script to run multiple SLUG trials using MPI. Basic usage is:

mpirun -np 1 slug_mpi.py PARAMFILE --nproc NPROC --batchsize BATCHSIZE

where PARAMFILE is the name of the parameter file for the run, NPROC
is the number of MPI processes to use, and BATCHSIZE is the number of
slug trials per MPI process. Note that, even when using multiple MPI
processes (i.e., when NPROC is > 1), it is important to call mpirun
with the argument -np 1. This is because mpirun should launch only 1
instance of the master script, and the master script will then launch
the appropriate number of MPI tasks on its own.
"""

import argparse
import copy
from mpi4py import MPI
import numpy as np
import os
import os.path as osp
import subprocess
import threading
import sys
import warnings
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

    
# Step 1: grab command line arguments and set up directory pointers
parser = argparse. \
         ArgumentParser(
             description="Wrapper script to run slug using mpi parallelism")
parser.add_argument('paramfile', help="name of slug parameter file")
parser.add_argument('-n', '--nproc', default=None, type=int,
                    help="number of slug processes (default: "+
                    "number of cores)")
parser.add_argument('-b', '--batchsize', default=None, type=int,
                    help="number of trials per slug process "+
                    "(default: ntrial/nproc)")
parser.add_argument('-nc', '--noconsolidate', action='store_true',
                    default=False, help="leave outputs in separate "
                    "files (default action: consolidate into a single "
                    "file)")
parser.add_argument('-o', '--outdir', default=None, type=str,
                    help="directory into which to write output; "
                    "defaults to current directory")
parser.add_argument('-r', '--restart', default=False,
                    action='store_true',
                    help="indicates that this is a restart of a "
                    "previous run for which some trials have already "
                    "been completed")
parser.add_argument('-v', '--verbose', action='store_true',
                    default=False, help="produce verbose output")
def return_parser():
    # Function to return the parser -- used by sphinx autodoc
    return parser
args = parser.parse_args()
scrdir = osp.dirname(osp.realpath(__file__))
cwd = os.getcwd()
if args.outdir is None:
    outdir = cwd
else:
    outdir = args.outdir
tmpdir = osp.join(outdir, 'tmp')    

# Step 2: suck the slug paramter file into memory
try:
    fp = open(args.paramfile, 'r')
except IOError:
    fp = None
    if not osp.isabs(args.paramfile) and 'SLUG_DIR' in os.environ:
        fname = osp.join(os.environ['SLUG_DIR'], args.paramfile)
        try:
            fp = open(fname, 'r')
        except IOError:
            pass
if fp is None:
    raise IOError("slug error: unable to open file "+args.paramfile)
pfile = fp.read().split('\n')
fp.close()


# Step 3: parse the file to grab the data we need out of it
ntrials = 1
ntrials_line = -1
model_name = 'SLUG_DEF'
model_name_line = -1
out_dir = cwd
out_dir_line = -1
output_mode = 'ascii'
verbosity = 0
sim_type = 'galaxy'
rng_offset_line = -1
for i, line in enumerate(pfile):
    linesplit = line.split()
    if len(linesplit) == 0:
        continue
    if linesplit[0].lower() == 'n_trials':
        try:
            ntrials = int(linesplit[1])
            ntrials_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'model_name':
        try:
            model_name = linesplit[1]
            model_name_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'out_dir':
        try:
            out_dir = linesplit[1]
            if not osp.isabs(out_dir):
                out_dir = osp.join(cwd, out_dir)
            out_dir_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'output_mode':
        try:
            output_mode = linesplit[1]
            output_mode = output_mode.lower()
            if output_mode == 'ascii':
                extension = '.txt'
            elif output_mode == 'binary':
                extension = '.bin'
            elif output_mode == 'fits':
                extension = '.fits'
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'rng_offset':
        try:
            rng_offset_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'verbosity':
        try:
            verbosity = int(linesplit[1])
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'sim_type':
        try:
            sim_type = linesplit[1]
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)

        
# Step 4: if this is a restart, search for already-existing
# output files and see how many completed trials they contain
completed_trials = 0
task_ctr = 0
if args.restart:
    while True:
        # Look for a summary file
        fname = osp.join(out_dir,
                         model_name+"_{:05d}_summary.txt".
                         format(task_ctr))
        if not osp.isfile(fname):
            break

        # Read number of trials from file
        fp = open(fname, 'r')
        for line in fp:
            linesplit = line.split()
            if len(linesplit) == 0:
                continue
            if linesplit[0].lower() == 'n_trials':
                completed_trials += int(linesplit[1])
                break
        fp.close()
        task_ctr += 1

# Step 5: set number of processors and batch size
if args.nproc is None:
    # Number of processes not specified, so try to figure it out from
    # environment variables
    if 'PBS_NCPUS' in os.environ:
        # PBS / torque
        nproc = int(os.environ['PBS_NCPUS']) - 1
    elif 'SLURM_JOB_NUM_NODES' in os.environ:
        # SLURM
        nproc = int(os.environ['SLURM_JOB_NUM_NODES']) * \
                int(os.environ['SLURM_CPUS_ON_NODE']) - 1
    elif 'LOADL_PROCESSOR_LIST' in os.environ:
        # LoadLeveler
        nproc = len(os.environ['LOADL_PROCESSOR_LIST'].split()) - 1
    else:
        raise RuntimeError('unable to determine CPU count from '
                           'environment variables; set manually '
                           'using --nprocs')
else:
    nproc = args.nproc
if args.batchsize is None:
    batchsize = int(np.ceil(float(ntrials-completed_trials)/nproc))
    if batchsize == 0:
        batchsize = 1
else:
    batchsize = args.batchsize
tot_jobs = task_ctr + \
           int(np.ceil(float(ntrials-completed_trials) / batchsize))

# Step 6: make temporary directory into which outputs can be written
try:
    os.mkdir(tmpdir)
except OSError: pass           # Probably failed because dir exists


# Step 7: main loop; we loop over tasks to be done, launching slug
# instances as needed; we need holders for the MPI communicators with
# children,
comms = [None]*nproc  # Holder for communicators with childern
req = [None]*nproc    # Holder for MPI requests
jobnum = np.zeros(nproc, dtype='int') # Job numbers
recv_buf = np.zeros(1, dtype='int')   # Buffer for disconnect messages
assigned_trials = np.zeros(nproc, dtype='int')
while completed_trials < ntrials:

    # Step 7a: check for non-busy tasks, and assign work
    for p in range(nproc):
        if assigned_trials[p] == 0:

            # Non-busy task, so assign work

            # Increment assigned trials list and task counter
            assigned_trials[p] \
                = min(batchsize,
                      ntrials - completed_trials -
                      np.sum(assigned_trials))
            if assigned_trials[p] == 0:
                continue
            task_ctr = task_ctr + 1
            jobnum[p] = task_ctr

            # Create a new slug parameter file with different number of
            # trials, output name, output directory, and rng offset
            pfile_tmp = copy.deepcopy(pfile)
            if ntrials_line != -1:
                pfile_tmp[ntrials_line] \
                    = "n_trials   {:d}".format(ntrials)
            else:
                pfile_tmp.append("n_trials   {:d}".format(ntrials))
            if model_name_line != -1:
                pfile_tmp[model_name_line] \
                    = "model_name   "+model_name+"_{:05d}".format(jobnum)
            else:
                pfile_tmp.append("model_name   "+model_name+
                                 "_{:05d}".format(jobnum))
            if out_dir_line != -1:
                pfile_tmp[out_dir_line] = "out_dir   "+tmpdir
            else:
                pfile_tmp.append("out_dir   "+tmpdir)
            if rng_offset_line != -1:
                pfile_tmp[rng_offset_line] \
                    = "rng_offset   " + str(10000*jobnum)
            else:
                pfile_tmp.append("rng_offset   " + str(10000*jobnum))
            pfile_name = osp.join(tmpdir,
                                  'slug_par_{:05d}'.format(jobnum))
            fp = open(pfile_name, 'w')
            for line in pfile_tmp:
                fp.write(line+'\n')
            fp.close()

            # Print status if verbose
            if args.verbose:
                print(("mpi rank {:d}: launching slug run " +
                       "{:d} / {:d} with {:d} trials").
                      format(rank, jobnum+1, tot_jobs, ntrials))

            # Launch a slug process
            slugcmd = osp.join(scrdir, '..', 'bin', 'slug')
            slugargs = [pfile_name, '--MPI_worker']
            comms[p] = MPI.COMM_SELF.Spawn(slugcmd, args=slugargs,
                                           maxprocs=1)
            
            # Register a non-blocking receive with this process;
            # it will receive the all done message
            if assigned_trials[p] > 0:
                req[p] = comm.Irecv(recv_buf, source=p, tag=2)

    # Step 7b: check for completion of running jobs
    for p in range(nproc):
        if req[p] is not None:
            if req[p].Test():

                # A slug job has completed

                # Move output files from temporary directory to main
                # directory
                basename = model_name+"_{:05d}".format(jobnum[p])
                endings = ["summary.txt",
                           "integrated_prop"+extension,
                           "integrated_phot"+extension,
                           "integrated_spec"+extension,
                           "integrated_yield"+extension,
                           "cluster_prop"+extension,
                           "cluster_phot"+extension,
                           "cluster_spec"+extension,
                           "cluster_yield"+extension]
                for e in endings:
                    fname = osp.join(tmpdir, osp.basename(basename+"_"+e))
                    if osp.isfile(fname):
                        os.rename(fname, osp.join(out_dir,
                                                  osp.basename(fname)))

                # Delete parameter file
                os.remove(pfile_name)

                # Update task information
                completed_trials += assigned_trials[p]
                assigned_trials[p] = 0
                req[p] = None
                jobnum[p] = -1
                comms[p].Disconnect()

                # Print status
                if args.verbose:
                    print("... {:d} / {:d} trials done".format(
                        completed_trials, ntrials))
        
# Step 8: consolidate if requested
if args.noconsolidate == False:

    combined_name = osp.join(out_dir, model_name)
    if verbosity > 0:
        print("Consolidating outputs to " + combined_name)

    # Step 8a: summary files: change the model name, output
    # directory, and number of trials in the first output, and
    # delete all the rest
    fp = open(osp.join(out_dir, model_name+'_00000_summary.txt'), 'r')
    fpout = open(osp.join(out_dir, model_name+'_summary.txt'), 'w')
    for line in fp:
        linesplit = line.split()
        if linesplit[0] == 'model_name':
            fpout.write("model_name           "+model_name+"\n")
        elif linesplit[0] == 'out_dir':
            fpout.write("out_dir              "+out_dir+"\n")
        elif linesplit[0] == 'n_trials':
            fpout.write("n_trials             {:d}\n".format(ntrials))
        else:
            fpout.write(line)
    fp.close()
    fpout.close()

    # Step 8b: integrated files: read data from all files, combine,
    # then write back out
    if sim_type != 'cluster':
        data = []
        nointegrated = False
        for i in range(tot_jobs):
            f = osp.join(out_dir,
                         model_name+"_{:05d}".format(i))
            if verbosity > 1:
                print("Reading integrated data from "+f+"...")
            try:
                data.append(read_integrated(f, fmt=output_mode,
                                            nofilterdata=True))
            except IOError:
                nointegrated = True
                break
        if not nointegrated:
            combined_data = combine_integrated(data)
            write_integrated(combined_data, combined_name,
                             fmt=output_mode)

    # Step 8c: cluster files; same as integrated files
    data = []
    nocluster = False
    for i in range(tot_jobs):
        f = osp.join(out_dir,
                     model_name+"_{:05d}".format(i))
        if verbosity > 1:
            print("Reading cluster data from "+f+"...")
        try:
            data.append(read_cluster(f, fmt=output_mode,
                                     nofilterdata=True))
        except IOError:
            nocluster = True
            break
    if not nocluster:
        combined_data = combine_cluster(data)
        write_cluster(combined_data, combined_name, fmt=output_mode)

    # Step 8d: clean up remaining temporary files
    if verbosity > 0:
        print("Cleaning up temporary files")
    endings = ["summary.txt",
               "integrated_prop"+extension,
               "integrated_phot"+extension,
               "integrated_spec"+extension,
               "integrated_yield"+extension,
               "cluster_prop"+extension,
               "cluster_phot"+extension,
               "cluster_spec"+extension,
               "cluster_yield"+extension]
    for i in range(tot_jobs):
        for e in endings:
            f = osp.join(out_dir,
                         model_name+"_{:05d}".format(i)+"_"+e)
            try:
                if osp.isfile(f):
                    os.remove(f)
            except OSError:
                warnings.warn("unable to clean up temporary file "+f)

# Step 9: remove temporary directory
try:
    os.rmdir(tmpdir)
except OSError:
    warnings.warn("unable to clean up temporary directory "+
                  tmpdir)
   

        # Send a message to parent stating that we have completed work
        job_done = np.ones(1, dtype='int')
        comm.Send([job_done, MPI.INT], dest=0, tag=2)

    # Step 7w: disconnect when we have been told to do so
    comm.Disconnect()
