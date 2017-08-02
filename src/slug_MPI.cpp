
#include <cstdlib>
#include <stdexcept>
#include "slug_MPI.H"

using namespace std;

////////////////////////////////////////////////////////////////////////
// An abort routine
////////////////////////////////////////////////////////////////////////
[[noreturn]] void bailout(int exit_val) {
#ifdef ENABLE_MPI
  MPI_Abort(MPI_COMM_WORLD, exit_val);
#endif
  exit(exit_val);
}

#ifdef ENABLE_MPI

////////////////////////////////////////////////////////////////////////
// Routines for exchanging clusters via MPI
////////////////////////////////////////////////////////////////////////

// Default max chunk size = 4 MB; this is set by a tradeoff of
// efficiency vs. memory. We want to use big chunks so as to guarantee
// that we send as few messages as possible, but we need to allocate a
// minimum memory size of SLUG_MPI_CHUNK_SIZE on each process, so
// big chunk sizes can be memory expensive.
#define SLUG_MPI_CHUNK_SIZE 4194304

// Utility routine for creating a send buffer from a vector of
// clusters
static inline slug_cluster_buffer *
pack_slug_clusters(const vector<slug_cluster *> &clusters,
		   vector<size_t> &sizes, size_t &bufsize) {
  
  // Record size of each cluster and total buffer size
  sizes.resize(clusters.size());
  bufsize = 0;
  for (vector<size_t>::size_type i=0; i<clusters.size(); i++) {
    sizes[i] = clusters[i]->buffer_size();
    bufsize += sizes[i];
  }

  // Allocate and pack the buffer
  slug_cluster_buffer *buf = malloc(bufsize);
  size_t ptr = 0;
  for (vector<size_t>::size_type i=0; i<clusters.size(); i++) {
    clusters[i]->pack_buffer((char *) buf+ptr);
    ptr += sizes[i];
  }

  // Return the buffer
  return buf;
}

// Utility routine to unpack a vector of clusters
static inline vector<slug_cluster *>
unpack_slug_clusters(const vector<size_t>::size_type ncluster,
		     const size_t *sizes,
		     const slug_cluster_buffer *buf,
		     const slug_PDF *imf_, 
		     const slug_tracks *tracks_, 
		     const slug_specsyn *specsyn_,
		     const slug_filter_set *filters_,
		     const slug_extinction *extinct_,
		     const slug_nebular *nebular_,
		     const slug_yields *yields_,
		     const slug_line_list *lines_,
		     slug_ostreams &ostreams_,
		     const slug_PDF *clf_) {
  vector<slug_cluster *> clusters(ncluster);
  size_t ptr = 0;
  for (vector<size_t>::size_type i=0; i<ncluster; i++) {
    slug_cluster_buffer *bufptr = (slug_cluster_buffer *)
      ((char *) buf + ptr);
    clusters[i] = new slug_cluster(bufptr, imf_, tracks_, specsyn_,
				   filters_, extinct_, nebular_,
				   yields_, lines_, ostreams_, clf_);
    ptr += sizes[i];
  }
  return clusters;
}

// Utility routine to pack a vector of clusters for chunked
// communication; this differs from pack_slug_clusters in that the
// buffer also contains information on the number of clusters and
// their sizes, packed into a single message. The buffer is broken up
// into chunks of a specified maximum size, and the routine returns
// the size of each chunk. The first chunk starts with
//
// total number of chunks [sizeof(size_t) bytes]
//
// number of clusters in each chunk [nchunk * sizeof(size_t) bytes]
//
// This is followed by:
//
// size of each cluster [ncluster * sizeof(size_t) bytes]
//
// data for each cluster
//
// All chunks but the first only contain the size of each cluster and
// the data for each cluster.
//
// For most applications the idea is that there will be only a single
// chunk, but for generality we cannot guarantee the data will all fit
// in one chunk, so we allow an arbitrary number
//
static slug_cluster_buffer *
pack_slug_clusters_chunk(const vector<slug_cluster *> &clusters,
			 vector<size_t> &chunk_sizes) {

  // Get size of each cluster to be packed; runtime error if any
  // individual cluster is too large to fit within maximum chunk size
  vector<size_t> sizes(clusters.size());
  for (vector<size_t>::size_type i=0; i<sizes.size(); i++) {
    sizes[i] = clusters[i]->buffer_size();
    if (sizes[i] > SLUG_MPI_CHUNK_SIZE)
      throw runtime_error("cluster memory size exceeds "
			  "SLUG_MPI_CHUNK_SIZE"
			  "in pack_slug_clusters_chunk; try recompiling"
			  "with a larger value for "
			  "SLUG_MPI_CHUNK_SIZE");
  }
      
  // Loop through the clusters, figuring out how to chunk them; note
  // that we traverse the cluster list backwards, because we don't
  // know in advance how much space we're going to need in the first
  // chunk to store the list of individual chunk sizes. Starting at
  // the back of the list guarantees that we minimise the odds that
  // we'll need to add an extra chunk to accommodate the header
  // information.
  vector<size_t> ncluster;
  chunk_sizes.resize(1);
  chunk_sizes[0] = 0;
  ncluster.push_back(0);
  
  // Loop through the clusters, figuring out how to chunk them
  for (long i=sizes.size()-1; i>=0; i--) {

    // See if we have room for this cluster in the current chunk
    if (chunk_sizes.front() + sizeof(size_t) + sizes[i]
	<= SLUG_MPI_CHUNK_SIZE) {
      // Yes, there is room, so increment size and number of clusters
      // in current chunk
      chunk_sizes.front() += sizes[i] + sizeof(size_t);
      ncluster.front()++;
    } else {
      // Current chunk is full; add a new chunk holder
      chunk_sizes.insert(chunk_sizes.begin(), sizes[i] + sizeof(size_t));
      ncluster.insert(ncluster.begin(), 1);
    }
  }

  // Now add size of header data to first chunk; if that makes it too
  // big, split the first chunk in two.
  chunk_sizes[0] += sizeof(size_t) * (chunk_sizes.size() + 1);
  if (chunk_sizes[0] > SLUG_MPI_CHUNK_SIZE) {
    // Create a new chunk in second position
    chunk_sizes.insert(chunk_sizes.begin()+1, 0);
    ncluster.insert(ncluster.begin()+1, 0);
    // Move clusters from the first chunk into the newly-created
    // second one until the size of the first chunk is small enough
    long i=ncluster[0];
    while (chunk_sizes[0] > SLUG_MPI_CHUNK_SIZE) {
      ncluster[1]++;
      ncluster[0]--;
      chunk_sizes[1] += sizes[i] + sizeof(size_t);
      chunk_sizes[0] -= sizes[i] + sizeof(size_t);
      i--;
      if (i < 0)
	throw runtime_error("insufficient memory chunk size in "
			    "pack_slug_clusters_chunk; try recompiling "
			    "with a larger value for "
			    "SLUG_MPI_CHUNK_SIZE");
    }
  }	

  // We now know how the data should be chunked; allocate a buffer
  // large enough to hold all the chunks
  size_t bufsize = 0;
  for (vector<size_t>::size_type i=0; i<chunk_sizes.size(); i++)
    bufsize += chunk_sizes[i];
  slug_cluster_buffer *buf = malloc(bufsize);

  // Store the total number of chunks and size of each chunk
  *((size_t *) buf) = (size_t) chunk_sizes.size();
  char *ptr = (char *) buf + sizeof(size_t);
  for (vector<size_t>::size_type i=0; i<ncluster.size(); i++) {
    *((size_t *) ptr) = ncluster[i];
    ptr += sizeof(size_t);
  }
  
  // Now fill in the sizes of each cluster and the cluster data
  vector<size_t>::size_type cluster_ptr = 0;
  for (vector<size_t>::size_type i=0; i<chunk_sizes.size(); i++) {

    // Fill in the sizes of the clusters in this chunk
    for (vector<size_t>::size_type j=0; j<ncluster[i]; j++) {
      *((size_t *) ptr) = sizes[cluster_ptr+j];
      ptr += sizeof(size_t);
    }

    // Write the cluster data for this chunk
    for (vector<size_t>::size_type j=0; j<ncluster[i]; j++) {
      clusters[cluster_ptr+j]->pack_buffer(ptr);
      ptr += sizes[cluster_ptr+j];
    }

    // Advance the cluster pointer
    cluster_ptr += ncluster[i];
  }

  // Return
  return buf;
}


// Blocking send of a single cluster
void MPI_send_slug_cluster(const slug_cluster &cluster, int dest, int tag,
			   MPI_Comm comm) {
  
  // Get size of buffer
  size_t bufsize = cluster.buffer_size();

  // Send size of buffer; note that this is a bit tricky, because
  // there isn't a single MPI integer type that we can guarantee will
  // correspond to size_t; we therefore send this as MPI_BYTE, and
  // interpret on the other end
  MPI_Send(&bufsize, sizeof(bufsize), MPI_BYTE, dest, tag, comm);

  // Get the serialized version of the input object
  slug_cluster_buffer *buf = cluster.make_buffer();

  // Do send and wait for completion
  MPI_Send(buf, bufsize, MPI_BYTE, dest, tag, comm);

  // Free buffer
  cluster.free_buffer(buf);
}

// Blocking receive of a single cluster
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
		      const slug_PDF *clf_) {

  // Receive the buffer size we'll need
  size_t bufsize;
  MPI_Recv(&bufsize, sizeof(bufsize), MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Allocate a buffer to hold the data
  slug_cluster_buffer *buf = malloc(bufsize);

  // Receive the cluster buffer
  MPI_Recv(buf, bufsize, MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Build a new slug_cluster object from the buffer
  slug_cluster *cluster
    = new slug_cluster(buf, imf_, tracks_, specsyn_, filters_,
		       extinct_, nebular_, yields_, lines_, ostreams_, clf_);

  // Free the buffer
  free(buf);

  // Return
  return cluster;
}

// Blocking send of a vector of clusters
void MPI_send_slug_cluster_vec(const vector<slug_cluster *> &clusters,
			       int dest, int tag, MPI_Comm comm) {

  // Create send buffer and metadata
  vector<size_t> sizes;
  size_t bufsize;
  slug_cluster_buffer *buf =
    pack_slug_clusters(clusters, sizes, bufsize);

  // First send number of objects to be sent
  vector<int>::size_type ncluster = sizes.size();
  MPI_Send(&ncluster, sizeof(ncluster), MPI_BYTE, dest, tag, comm);

  // Now send vector of sizes of each object
  MPI_Send(sizes.data(), sizeof(size_t)*ncluster, MPI_BYTE,
	   dest, tag, comm);

  // Finally send the data
  MPI_Send(buf, bufsize, MPI_BYTE, dest, tag, comm);

  // Free buffer
  free(buf);
}

// Blocking receive of a vector of clusters
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
			  const slug_PDF *clf_) {

  // First receive the number of objects to be received
  vector<size_t>::size_type ncluster;
  MPI_Recv(&ncluster, sizeof(ncluster), MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Now receive the sizes of each object
  size_t *sizes = (size_t *) calloc(ncluster, sizeof(size_t));
  MPI_Recv(sizes, sizeof(size_t)*ncluster, MPI_BYTE, source, tag,
	   comm, MPI_STATUS_IGNORE);

  // Allocate memory to receive the cluster data
  size_t bufsize = 0;
  for (vector<int>::size_type i=0; i<ncluster; i++) bufsize += sizes[i];
  slug_cluster_buffer *buf = (slug_cluster_buffer *) malloc(bufsize);

  // Receive data
  MPI_Recv(buf, bufsize, MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Unpack data
  vector<slug_cluster *> clusters =
    unpack_slug_clusters(ncluster, sizes, buf, imf_, tracks_, specsyn_,
			 filters_, extinct_, nebular_, yields_,lines_,
			 ostreams_, clf_);

  // Free buffers
  free(buf);
  free(sizes);

  // Return
  return clusters;
}

// Broadcasting of a single cluster
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
		       const slug_PDF *clf_) {

  // What we do here depends on if we are the root processor
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if (myrank == root) {

    // I am the transmitting processor, so get the size of the buffer
    // and build the buffer
    size_t bufsize = cluster->buffer_size();
    slug_cluster_buffer *buf = cluster->make_buffer();

    // Send the size of the buffer
    MPI_Bcast(&bufsize, sizeof(bufsize), MPI_BYTE, root, comm);

    // Send the contents of the buffer
    MPI_Bcast(buf, bufsize, MPI_BYTE, root, comm);

    // Free buffer
    free(buf);

    // Return the input pointer as output
    return cluster;

  } else {

    // I am receiving

    // Receive the size of the buffer
    size_t bufsize;
    MPI_Bcast(&bufsize, sizeof(bufsize), MPI_BYTE, root, comm);

    // Allocate memory to hold the buffer, then receive it
    void *buf = malloc(bufsize);
    MPI_Bcast(buf, bufsize, MPI_BYTE, root, comm);

    // Construct the new slug_cluster object
    slug_cluster *new_cluster =
      new slug_cluster(buf, imf_, tracks_, specsyn_, filters_,
		       extinct_, nebular_, yields_,lines_, ostreams_, clf_);

    // Free the buffer
    free(buf);

    // Return the new object
    return new_cluster;
  }
}

// Broadcasting a vector of clusters
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
			   const slug_PDF *clf_) {

  // What we do here depends on if we are the root processor
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  if (myrank == root) {

    // I am the sending processor, so build the send buffer and its
    // metadata
    vector<size_t> sizes;
    size_t bufsize;
    slug_cluster_buffer *buf =
      pack_slug_clusters(clusters, sizes, bufsize);
    vector<size_t>::size_type ncluster = sizes.size();

    // Broadcast how many clusters we'll be sending
    MPI_Bcast(&ncluster, sizeof(ncluster), MPI_BYTE, root, comm);

    // Now broadcast the sizes of each cluster
    MPI_Bcast(sizes.data(), ncluster*sizeof(size_t), MPI_BYTE, root, comm);

    // Now broadcast the data
    MPI_Bcast(buf, bufsize, MPI_BYTE, root, comm);

    // Free buffer
    free(buf);

    // Return the array of pointers we were originally passed
    return clusters;

  } else {

    // I am a receiving processor

    // Receive the number of clusters
    vector<size_t>::size_type ncluster = 0;
    MPI_Bcast(&ncluster, sizeof(ncluster), MPI_BYTE, root, comm);

    // Receive the sizes of each individual cluster
    size_t *sizes = (size_t *) calloc(ncluster, sizeof(size_t));
    MPI_Bcast(sizes, ncluster*sizeof(size_t), MPI_BYTE, root, comm);

    // Allocate memory to receive the cluster data
    size_t bufsize = 0;
    for (vector<int>::size_type i=0; i<ncluster; i++) bufsize += sizes[i];
    slug_cluster_buffer *buf = (slug_cluster_buffer *) malloc(bufsize);

    // Receive the data
    MPI_Bcast(buf, bufsize, MPI_BYTE, root, comm);

    // Unpack data
    vector<slug_cluster *> clusters =
      unpack_slug_clusters(ncluster, sizes, buf, imf_, tracks_, specsyn_,
			   filters_, extinct_, nebular_, yields_, lines_,
			   ostreams_, clf_);

    // Free buffers
    free(buf);
    free(sizes);

    // Return
    return clusters;
  }
}

////////////////////////////////////////////////////////////////////////
// Exchanging clusters
////////////////////////////////////////////////////////////////////////

void
MPI_exchange_slug_cluster(const vector<slug_cluster *> &clusters,
			  const vector<int> &destinations,
			  vector<slug_cluster *> &received_clusters,
			  vector<int> &sources,
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
			  const slug_PDF *clf_) {

  // Get MPI information
  int nrank, myrank;
  MPI_Comm_size(comm, &nrank);
  MPI_Comm_rank(comm, &myrank);

  // Prepare storage to hold the data we'll be sending and receiving
  vector<vector<slug_cluster *> > send_clusters(nrank),
    recv_clusters(nrank);
  vector<slug_cluster_buffer *> send_buf(nrank);
  vector<vector<slug_cluster_buffer *> > recv_buf(nrank);
  vector<vector<size_t>::size_type> recv_nchunk(nrank);
  vector<vector<vector<size_t>::size_type> >  recv_ncluster(nrank);

  // Sort clusters by destination; for sends to myself, just copy the
  // pointers directly to the output list
  for (vector<int>::size_type i=0; i<destinations.size(); i++) {
    if (destinations[i] < 0) continue;
    else if (destinations[i] == myrank) {
      recv_clusters[i].push_back(clusters[i]);
    } else {
      send_clusters[destinations[i]].push_back(clusters[i]);
    }
  }

  // Record status of communication at start
  vector<bool> send_done(nrank), recv_done(nrank);
  vector<vector<MPI_Request> > send_req(nrank), recv_req(nrank);
  for (int i=0; i<nrank; i++) {
    if (i != myrank)
      send_done[i] = recv_done[i] = false;
    else
      send_done[i] = recv_done[i] = true;
  }

  // Enter main loop
  bool comm_done = false;
  while (!comm_done) {

    // Set communication done to true; will be changed below if not
    // done
    comm_done = true;

    // Loop over ranks
    for (int i=0; i < nrank; i++) {

      // Skip send-receive to myself
      if (i == myrank) continue;

      // Send block, executed if the send to this processor isn't done
      if (!send_done[i]) {

	// Send not done, so communication not done
	comm_done = false;

	// Have we already registered the send to this target?
	if (send_req[i].size() == 0) {

	  // We have not sent to this target yet, so make a buffer and
	  // register a non-blocking send of each chunk
	  vector<size_t> chunk_sizes;
	  send_buf[i] =
	    pack_slug_clusters_chunk(send_clusters[i], chunk_sizes);
	  send_req[i].resize(chunk_sizes.size());
	  char *ptr = (char *) send_buf[i];
	  for (vector<size_t>::size_type j=0; j<chunk_sizes.size(); j++) {
	    MPI_Isend(ptr, chunk_sizes[j], MPI_BYTE, i, j, comm,
		      &(send_req[i][j]));
	    ptr += chunk_sizes[j];
	  }

	} else {

	  // We have sent to this target; check if all of our requests
	  // have completed, and, if so, free the buffer and record that
	  // we are done
	  bool all_done = true;
	  for (vector<MPI_Request>::size_type j=0; j<send_req[i].size();
	       j++) {
	    
	    // Check pending incomplete requests
	    if (send_req[i][j] != MPI_REQUEST_NULL) {
	      int flag;
	      MPI_Status stat;
	      MPI_Test(&(send_req[i][j]), &flag, &stat);
	      all_done = all_done && flag;
	    }
	  }

	  // If all requests are done, free buffer and record that this
	  // exchange is done
	  if (all_done) {
	    free(send_buf[i]);
	    send_done[i] = true;
	  }
	  
	} // End testing for send completion

      } // End of send block

      // Receive block, executed if receive from this processor isn't
      // done
      if (!recv_done[i]) {

	// Receive not done, so all communication note done
	comm_done = false;

	// Have we already registered the send to this target?
	if (recv_req[i].size() == 0) {

	  // No receive request registered yet, so register a request
	  // to receive one chunk
	  recv_buf[i].push_back(malloc(SLUG_MPI_CHUNK_SIZE));
	  recv_req[i].push_back(MPI_REQUEST_NULL);
	  MPI_Irecv(recv_buf[i][0], SLUG_MPI_CHUNK_SIZE, MPI_BYTE,
		    i, 0, comm,
		    recv_req[i].data());
	  recv_nchunk[i] = 1;

	} else {

	  // We have registered at least one receive request; check
	  // for completion of outstanding receive requests
	  bool all_done = true;
	  for (vector<MPI_Request>::size_type j=0; j<recv_req[i].size();
	       j++) {
	    
	    // Check pending incomplete requests
	    if (recv_req[i][j] != MPI_REQUEST_NULL) {
	      int flag;
	      MPI_Status stat;
	      MPI_Test(&(recv_req[i][j]), &flag, &stat);
	      all_done = all_done && flag;

	      // If still pending, go on to next request
	      if (!flag) continue;

	      // If we just received the first chunk, see how many
	      // chunks in total there will be, and get the number of
	      // clusters in each; set the cluster array to the
	      // correct size
	      char *bufptr = (char *) recv_buf[i][j];
	      if (j==0) {
		recv_nchunk[i] = *((size_t *) bufptr);
		recv_ncluster[i].resize(recv_nchunk[i]);
		bufptr += sizeof(size_t);
		vector<int>::size_type ncluster_tot = 0;
		for (vector<int>::size_type k=0; k<recv_nchunk[i]; k++) {
		  recv_ncluster[i][k] = *((size_t *) bufptr);
		  bufptr += sizeof(size_t);
		  ncluster_tot += recv_ncluster[i][k];
		}
		recv_clusters[i].resize(ncluster_tot);
		if (ncluster_tot == 0) break;
	      }

	      // Unpack this chunk into the output array
	      size_t *sizes = (size_t *) bufptr;
	      bufptr += recv_ncluster[i][j] * sizeof(size_t);
	      vector<int>::size_type start = 0, stop;
	      for (vector<int>::size_type k=0; k<j; k++)
		start += recv_ncluster[i][k];
	      stop = start + recv_ncluster[i][j];
	      for (vector<int>::size_type k=start; k<stop; k++) {
		recv_clusters[i][k] =
		  new slug_cluster((slug_cluster_buffer *) bufptr,
				   imf_, tracks_, specsyn_,
				   filters_, extinct_, nebular_,
				   yields_, lines_, ostreams_, clf_);
		bufptr += sizes[k-start];
	      }

	      // Free the buffer
	      free(recv_buf[i][j]);
	    }
	  }

	  // Do we need to register more receive requests because the
	  // first chunk we received says there are more chunks? If
	  // so, do that now.
	  if (recv_nchunk[i] > recv_req[i].size()) {
	    recv_buf[i].resize(recv_nchunk[i]);
	    recv_req[i].resize(recv_nchunk[i]);
	    for (vector<MPI_Request>::size_type j=1; j<recv_nchunk[i];
		 j++) {
	      recv_buf[i][j] = malloc(SLUG_MPI_CHUNK_SIZE);
	      MPI_Irecv(recv_buf[i][j], SLUG_MPI_CHUNK_SIZE, MPI_BYTE,
			i, j, comm,
			recv_req[i].data() + j);
	    }
	    all_done = false;
	  }

	  // If we're all done, record that
	  if (all_done) recv_done[i] = true;

	} // End testing for receive completion

      } // End receive block
      
    } // End loop over ranks

  } // End main loop

  // Pack the clusters we've received into the final output array, and
  // record the source for each cluster we're returning
  vector<int>::size_type ncluster_tot = 0, ptr = 0;
  for (int i=0; i<nrank; i++) ncluster_tot += recv_clusters[i].size();
  received_clusters.resize(ncluster_tot);
  sources.resize(ncluster_tot);
  for (int i=0; i<nrank; i++) {
    for (vector<int>::size_type j=0; j<recv_clusters[i].size(); j++) {
      received_clusters[ptr] = recv_clusters[i][j];
      sources[ptr] = i;
      ptr += 1;
    }
  }
}


#endif
// ENABLE_MPI
