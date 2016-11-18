
#ifdef ENABLE_MPI

#include "slug_MPI.H"

using namespace std;

// Utility routine for creating a send buffer from a vector of
// clusters
static inline slug_cluster_buffer *
pack_slug_clusters(const vector<slug_cluster *> &clusters,
		   std::vector<size_t> &sizes, size_t &bufsize) {
  
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
		     const slug_PDF *clf_) {
  vector<slug_cluster *> clusters(ncluster);
  size_t ptr = 0;
  for (vector<size_t>::size_type i=0; i<ncluster; i++) {
    slug_cluster_buffer *bufptr = (slug_cluster_buffer *)
      ((char *) buf + ptr);
    clusters[i] = new slug_cluster(bufptr, imf_, tracks_, specsyn_,
				   filters_, extinct_, nebular_,
				   yields_, clf_);
    ptr += sizes[i];
  }
  return clusters;
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
		       extinct_, nebular_, yields_, clf_);

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
			 filters_, extinct_, nebular_, yields_, clf_);

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
		       extinct_, nebular_, yields_, clf_);

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
			   filters_, extinct_, nebular_, yields_, clf_);

    // Free buffers
    free(buf);
    free(sizes);

    // Return
    return clusters;
  }
}

#endif
// ENABLE_MPI
