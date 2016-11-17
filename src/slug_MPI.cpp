
#ifdef ENABLE_MPI

#include "slug_MPI.H"

using namespace std;

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

// Blocking send of a vector of clusters
void MPI_send_slug_cluster_vec(const vector<slug_cluster *> &clusters,
			       int dest, int tag, MPI_Comm comm) {

  // Record number of clusters and buffer sizes
  vector<size_t> bufsize(clusters.size());
  size_t bufsize_sum = 0;
  for (vector<int>::size_type i=0; i<clusters.size(); i++) {
    bufsize[i] = clusters[i]->buffer_size();
    bufsize_sum += bufsize[i];
  }

  // Allocate and pack the buffer
  slug_cluster_buffer *buf = malloc(bufsize_sum);
  size_t ptr = 0;
  for (vector<int>::size_type i=0; i<clusters.size(); i++) {
    clusters[i]->pack_buffer((char *) buf+ptr);
    ptr += bufsize[i];
  }

  // First send number of objects to be sent
  vector<int>::size_type ncluster = bufsize.size();
  MPI_Send(&ncluster, sizeof(ncluster), MPI_BYTE, dest, tag, comm);

  // Now send vector of sizes of each object
  MPI_Send(bufsize.data(), sizeof(size_t)*bufsize.size(), MPI_BYTE,
	   dest, tag, comm);

  // Finally send the data
  MPI_Send(buf, bufsize_sum, MPI_BYTE, dest, tag, comm);

  // Free buffer
  free(buf);
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
  vector<int>::size_type ncluster;
  MPI_Recv(&ncluster, sizeof(ncluster), MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Now receive the sizes of each object
  size_t *sizes = (size_t *) calloc(ncluster, sizeof(size_t));
  MPI_Recv(sizes, sizeof(size_t)*ncluster, MPI_BYTE, source, tag,
	   comm, MPI_STATUS_IGNORE);

  // Allocate memory to receive the cluster data
  size_t bufsize = 0;
  for (vector<int>::size_type i=0; i<ncluster; i++) bufsize += sizes[i];
  slug_cluster_buffer *buf = malloc(bufsize);

  // Receive data
  MPI_Recv(buf, bufsize, MPI_BYTE, source, tag, comm,
	   MPI_STATUS_IGNORE);

  // Construct slug cluster particles to return
  vector<slug_cluster *> clusters(ncluster);
  size_t ptr = 0;
  for (vector<int>::size_type i=0; i<ncluster; i++) {
    slug_cluster_buffer *bufptr = (slug_cluster_buffer *)
      ((char *) buf + ptr);
    clusters[i] = new slug_cluster(bufptr, imf_, tracks_, specsyn_,
				   filters_, extinct_, nebular_,
				   yields_, clf_);
    ptr += sizes[i];
  }

  // Free buffers
  free(buf);
  free(sizes);

  // Return
  return clusters;
}
#endif
// ENABLE_MPI
