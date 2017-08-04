/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

////////////////////////////////////////////////////////////////////////
//
// This module defines MPI bindings for slug. They make it possible to
// exchange slug objects via MPI. At present bindings are provided only
// to transmit slug_cluster objects.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_MPI_H_
#define _slug_MPI_H_

////////////////////////////////////////////////////////////////////////
// An abort routine
////////////////////////////////////////////////////////////////////////
[[noreturn]] void bailout(int exit_val);

#ifdef ENABLE_MPI

////////////////////////////////////////////////////////////////////////
// If we are trying to use an old version of MPI that doesn't allow
// all functionality, emit a compiler warning.
////////////////////////////////////////////////////////////////////////
#if (MPI_VERSION == 1) || (MPI_VERSION == 2)
#   ifndef _MPI_WARNING_PRINTED_
#      define _MPI_WARNING_PRINTED_
#      warning "Warning: you appear to be using an older MPI that does not support RMA functionality. SLUG will operate without problems as a library under this MPI version, but parallel functionality of the slug executable will be disabled."
#   endif
#endif



#include <iostream>
#include <vector>
#include "mpi.h"
#include "slug.H"
#include "slug_cluster.H"

////////////////////////////////////////////////////////////////////////
// Blocking send/receive of single clusters and vectors of clusters
////////////////////////////////////////////////////////////////////////

// Sending processor passes in a cluster and a target to which it is
// sent
void MPI_send_slug_cluster(const slug_cluster &cluster, int dest, int tag,
			   MPI_Comm comm);

// Receiving processor specifies who the sender is, and gets back a
// pointer to a newly-constructed slug_cluster and matches the one
// that was sent
slug_cluster *
MPI_recv_slug_cluster(int source, int tag, MPI_Comm comm,
		      const slug_PDF *imf_, 
		      const slug_tracks *tracks_, 
		      const slug_specsyn *specsyn_,
		      const slug_filter_set *filters_,
		      const slug_extinction *extinct_,
		      const slug_nebular *nebular_,
		      const slug_yields *yields_,
		      const slug_line_list *lines_,
		      slug_ostreams &ostreams_,
		      const slug_PDF *clf_=NULL);

// Sender passes in a vector of pointers to slug_cluster, which are
// all sent to a destination processor
void MPI_send_slug_cluster_vec(const std::vector<slug_cluster *> &clusters,
			       int dest, int tag, MPI_Comm comm);

// Receiver specifies which processor to get data from; the returned
// value is a vector of pointers to slug_cluster; the pointers point
// to newly-constructed slug_clusters that are copies of the ones sent
std::vector<slug_cluster *>
MPI_recv_slug_cluster_vec(int source, int tag, MPI_Comm comm,
			  const slug_PDF *imf_, 
			  const slug_tracks *tracks_, 
			  const slug_specsyn *specsyn_,
			  const slug_filter_set *filters_,
			  const slug_extinction *extinct_,
			  const slug_nebular *nebular_,
			  const slug_yields *yields_,
			  const slug_line_list *lines_,
			  slug_ostreams &ostreams_,
			  const slug_PDF *clf_=NULL);


////////////////////////////////////////////////////////////////////////
// Broadcasting of single clusters and vectors of clusters
////////////////////////////////////////////////////////////////////////

// The broadcasting processor passes in a pointer to the cluster it
// wants to broadcast; receiving processors can pass in null. The
// return value is a pointer to a slug_cluster that is a copy of the
// one provided on the broadcasting processor. On the broadcasting
// processor, this pointer is a copy of the variable cluster, and on
// all other processors it is a pointer to a newly-constructed
// cluster that is a copy of the one that was sent.
slug_cluster *
MPI_bcast_slug_cluster(slug_cluster *cluster, int root,
		       MPI_Comm comm,
		       const slug_PDF *imf_, 
		       const slug_tracks *tracks_, 
		       const slug_specsyn *specsyn_,
		       const slug_filter_set *filters_,
		       const slug_extinction *extinct_,
		       const slug_nebular *nebular_,
		       const slug_yields *yields_,
		       const slug_line_list *lines_,
		       slug_ostreams &ostreams_,
		       const slug_PDF *clf_=NULL);

// The broadcasting processor passes in a vector of pointers to
// slug_clusters it wants to send; all other processors can pass an
// empty vector to clusters. The return value is a vector of pointers
// to slug_clusters. On the sending processor these pointers are just
// copies of the pointers in clusters, while on all other processors
// they are pointers to newly-constructed slug_cluster objects that
// are copies of the ones sent.
std::vector<slug_cluster *>
MPI_bcast_slug_cluster_vec(std::vector<slug_cluster *> &clusters,
			   int root, MPI_Comm comm,
			   const slug_PDF *imf_, 
			   const slug_tracks *tracks_, 
			   const slug_specsyn *specsyn_,
			   const slug_filter_set *filters_,
			   const slug_extinction *extinct_,
			   const slug_nebular *nebular_,
			   const slug_yields *yields_,
			   const slug_line_list *lines_,
			   slug_ostreams &ostreams_,
			   const slug_PDF *clf_=NULL);


////////////////////////////////////////////////////////////////////////
// Synchronisation of clusters. This is the most automated
// routine. The user invokes this routine on every processor,
// supplying a vector of the clusters to be sent and a destination for
// each cluster; destinations of -1 are skipped. The routine returns a
// vector of the clusters received from all other processors. This
// routine is blocking, but internal communications are non-blocking
// to maximise efficiency. In terms of function, the communication
// pattern in this routine is equivalent to MPI_Alltoallv. The routine
// preserves cluster order, in the sense that the clusters that are
// received are ordered by MPI rank (i.e., all clusters from processor
// 0 come first, those from processor 1 come next, etc.), and within a
// given MPI rank clusters are ordered as they were on the sending
// processor.
//
// Arguments:
//    clusters = vector of pointers to clusters to be sent
//    destinations = vector of MPI ranks giving desination for each
//                   cluster
//    received_clusters = vector of pointers to clusters received from
//                        other processors
//    sources = vector of MPI ranks listing the processor from which
//              each of the received clusters came
//    comm = MPI communicator to use
//    imf, tracks, specsyn, filters, extinct, nebular, yields, clf =
//       pointers to these objects, used to reconstruct clusters that
//       are received; not communicated
//
////////////////////////////////////////////////////////////////////////
void
MPI_exchange_slug_cluster(const std::vector<slug_cluster *> &clusters,
			  const std::vector<int> &destinations,
			  std::vector<slug_cluster *> &received_clusters,
			  std::vector<int> &sources,
			  MPI_Comm comm,
			  const slug_PDF *imf_, 
			  const slug_tracks *tracks_, 
			  const slug_specsyn *specsyn_,
			  const slug_filter_set *filters_,
			  const slug_extinction *extinct_,
			  const slug_nebular *nebular_,
			  const slug_yields *yields_,
			  const slug_line_list *lines_,
			  slug_ostreams &ostreams_,
			  const slug_PDF *clf_=NULL);

#endif
// ENABLE_MPI

#endif
// _slug_MPI_H_

