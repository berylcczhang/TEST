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

#ifdef ENABLE_MPI

#ifndef _slug_MPI_H_
#define _slug_MPI_H_

#include <vector>
#include "mpi.h"
#include "slug.H"
#include "slug_cluster.H"


////////////////////////////////////////////////////////////////////////
// Blocking send/receive of single clusters and vectors of clusters
////////////////////////////////////////////////////////////////////////

void MPI_send_slug_cluster(const slug_cluster &cluster, int dest, int tag,
			   MPI_Comm comm);
void MPI_send_slug_cluster_vec(const std::vector<slug_cluster *> &clusters,
			       int dest, int tag, MPI_Comm comm);

slug_cluster *
MPI_recv_slug_cluster(int source, int tag, MPI_Comm comm,
		      const slug_PDF *imf_, 
		      const slug_tracks *tracks_, 
		      const slug_specsyn *specsyn_,
		      const slug_filter_set *filters_,
		      const slug_extinction *extinct_,
		      const slug_nebular *nebular_,
		      const slug_yields *yields_,
		      const slug_PDF *clf_=NULL);

std::vector<slug_cluster *>
MPI_recv_slug_cluster_vec(int source, int tag, MPI_Comm comm,
			  const slug_PDF *imf_, 
			  const slug_tracks *tracks_, 
			  const slug_specsyn *specsyn_,
			  const slug_filter_set *filters_,
			  const slug_extinction *extinct_,
			  const slug_nebular *nebular_,
			  const slug_yields *yields_,
			  const slug_PDF *clf_=NULL);


#endif
// _slug_MPI_H_

#endif
// ENABLE_MPI

