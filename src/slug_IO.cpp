/*********************************************************************
Copyright (C) 2017 Robert da Silva, Michele Fumagalli, Mark Krumholz
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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
    typedef decltype(nullptr) nullptr_t;
}
#endif

#include "slug_IO.H"

////////////////////////////////////////////////////////////////////////
// slug_prefixbuf class
////////////////////////////////////////////////////////////////////////
int slug_prefixbuf::sync() { return this->sbuf->pubsync(); }

int slug_prefixbuf::overflow(int c) {
  if (c != std::char_traits<char>::eof()) {
    if (this->need_prefix
	&& !this->prefix.empty()
	&& this->prefix.size() != (unsigned int)
	this->sbuf->sputn(&this->prefix[0], this->prefix.size())) {
      return std::char_traits<char>::eof();
    }
    this->need_prefix = c == '\n';
  }
  return this->sbuf->sputc(c);
}


#ifdef ENABLE_MPI
////////////////////////////////////////////////////////////////////////
// MPI version of slug_oprefixstream
////////////////////////////////////////////////////////////////////////

// Constructor
slug_oprefixstream::
slug_oprefixstream(std::string const& prefix_, std::ostream& out,
		   MPI_Comm comm_, bool only_one_proc_, int tag_)
    : slug_prefixbuf(prefix_, out.rdbuf()), 
      std::ios(static_cast<std::streambuf*>(this)),
      std::ostream(static_cast<std::streambuf*>(this)),
      comm(comm_), tag(tag_), only_one_proc(only_one_proc_),
      target(out) {

  // Record rank
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);

  // Are we going to be writing output from only one proessor, or from
  // all?
  if (!only_one_proc && comm != MPI_COMM_NULL) {
    // If we are writing messages from all processors, we need to say
    // which processor is talking
    std::stringstream ss;
    ss << prefix << "process " << rank+1 << ": ";
    prefix = ss.str();
    // Pointer in outgoing buffer starts at 0
    outgoing_ptr = 0;
  }
}


// This method forces all communication to complete, and blocks until
// it does. The synchronization here is a bit tricky, because rank 0
// doesn't know when all the other processors are done. The way we
// handle this is that, when a processor other than 0 finds that all
// its pending requests have completed, it sends a null string to
// processor 0. Processor 0 continues looping, waiting for messages,
// until it receives null strings from all other processors.
void slug_oprefixstream::flush_communication() {

  // This method only does something if we are in MPI mode
  if (!only_one_proc && comm != MPI_COMM_NULL) {

    if (rank == 0) {
      
      // Rank 0: build structure to keep track of termination signals
      // we have received from other processors
      int nproc;
      MPI_Comm_size(comm, &nproc);
      if (nproc > 0) {
	std::vector<bool> flush_done(static_cast<unsigned int>(nproc-1));
	for (std::vector<bool>::size_type i=0; i<flush_done.size(); i++)
	  flush_done[i] = false;

	// Now loop to receive messages
	while (true) {

	  // See if we have received all termination messages; if so,
	  // exit the loop
	  bool all_done = true;
	  for (std::vector<bool>::size_type i=0; i<flush_done.size(); i++)
	    all_done = all_done && flush_done[i];
	  if (all_done) break;
	
	  // Do blocking receive until we get a message
	  MPI_Status stat;
	  char msg[SLUG_MAX_LINE_SIZE+1];
	  for (unsigned int i=0; i<SLUG_MAX_LINE_SIZE+1; i++) msg[i] = '\0';
	  MPI_Recv(msg, SLUG_MAX_LINE_SIZE, MPI_CHAR, MPI_ANY_SOURCE,
		   tag, comm, &stat);

	  // Check if message is a null string
	  if (msg[0] == '\0') {
	    // Yes, so mark that we have gotten the termination message
	    // from this processor
	    flush_done[stat.MPI_SOURCE-1] = true;
	  } else {
	    // No, so write received data to output stream
	    target << msg;
	  }
	}
      }

    } else {

      // All other ranks wait until all pending messages return
      if (req.size() > 0) {
	MPI_Status *stat = new MPI_Status[req.size()];
	MPI_Waitall(req.size(), req.data(), stat);
	delete [] stat;
      }

      // Erase our buffer and request list
      outgoing_buf.erase();
      outgoing_ptr = 0;
      req.resize(0);
      msg_size.resize(0);

      // Now send a null string to indicate completion
      char term_char = '\0';
      MPI_Send(&term_char, 1, MPI_CHAR, 0, tag, comm);

    }
  }
}  


// The overflow function: behavior depends on whether only_one_proc is
// set:
// if only_one_proc is set, then this function does nothing if called
// on any rank other than 0, and if called on rank 0, it acts like the
// serial verison
// if only_one_proc is not set, then on any rank other than 0, this
// function transfers characters to processor 0 for output; on
// processor 0, it grabs from the transfer buffer and writes it before
// writing its own next line of output
int slug_oprefixstream::overflow(int c) {

  if (comm == MPI_COMM_NULL) {

    // If we were given a null communicator, just act like the serial
    // version
    return slug_prefixbuf::overflow(c);

  } else if (only_one_proc) {
    
    // Case where processor 0 is printing only
    if (rank == 0) return slug_prefixbuf::overflow(c);
    else return 0;

  } else {
    
    // Case where we need to transfer data to rank 0, because all
    // processors want to write

    if (rank == 0) {
      // Action on rank 0

      // Process the character we got and store the return value
      int ret = slug_prefixbuf::overflow(c);

      // If the character we got is an endline, see if we have any
      // data transfers from other processors
      if (c == '\n') {

	// Loop until we have no pending messages
	while (true) {

	  // Check for pending messages
	  int message_avail;
	  MPI_Status stat;
	  MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &message_avail,
		     &stat);

	  // If no messages, stop looping
	  if (!message_avail) break;

	  // Message is available, so now receive it
	  char msg[SLUG_MAX_LINE_SIZE+1];
	  for (unsigned int i=0; i<SLUG_MAX_LINE_SIZE+1; i++) msg[i] = '\0';
	  MPI_Recv(msg, SLUG_MAX_LINE_SIZE, MPI_CHAR, stat.MPI_SOURCE,
		   tag, comm, &stat);

	  // Write received data to output stream
	  target << msg;
	}
      }
	
      // Return the result of writing our own line
      return ret;
      
    } else {

      // Action on all other ranks

      // First check all pending requests to see if they have
      // completed, in which case we can erase them from our outgoing
      // buffer
      while (req.size() > 0) {
	// Test to see if request has completed
	int req_complete;
	MPI_Status stat;
	MPI_Test(&req.front(), &req_complete, &stat);
	if (req_complete) {
	  // If yes, erase the corresponding elements from the outgoing
	  // buffer, and pop the request off the queue
	  outgoing_buf.erase(0, msg_size.front());
	  outgoing_ptr -= msg_size.front();
	  req.erase(req.begin());
	  msg_size.erase(msg_size.begin());
	} else {
	  // If no, stop checking, since requests should generally
	  // complete in the order they were posted
	  break;
	}
      }

      // Write prefix to our outgoing buffer if needed
      if (need_prefix) outgoing_buf += prefix;
      
      // Store the character we just got
      outgoing_buf += c;

      // Next action depends on whether the character we got is a
      // carriage return
      if (c != '\n') {

	// Character is not a carriage return, so add to outgoing buffer
	need_prefix = false;
	
      } else {

	// Character is a carriage return; mark this and add to our buffer
	need_prefix = true;

	// Now we do a little hack; pre-MPI 3.0, the C bindings for
	// MPI didn't declare the buffer argument to MPI_Isend as
	// const, even though it was, because you were expected to use
	// the C++ bindings, where it was declared as const. This was
	// changed in MPI 3.0, when the C++ bindings were deprecated,
	// but the net effect is that you have to be careful when
	// writing code that will function under either MPI
	// standard. The hack we use to solve this problem is that we
	// use the C bindings, since these work under any version of
	// MPI, but rather than using our string's buffer directly as
	// the argument to the send, which won't work under MPI 2.x
	// because the .c_str() method to string returns a const
	// buffer, we explicitly break const'ness by casting a pointer
	// to the buffer from const to non-const. This is safe,
	// because the function we then pass it to does in fact
	// respect const'ness, even if it doesn't say so in the
	// interface declaration under older versions of the MPI C
	// bindings... that was a lot of explanation for a one-line
	// statement...
	char *outgoing_buf_tmp 
	  = const_cast<char *>(outgoing_buf.c_str() + outgoing_ptr);
	
	// Now start a transfer
	unsigned int chars_to_transfer = outgoing_buf.size() - outgoing_ptr;
	req.push_back(MPI_REQUEST_NULL);
	msg_size.push_back(chars_to_transfer);
	MPI_Isend(outgoing_buf_tmp, chars_to_transfer, MPI_CHAR, 0, tag,
		  comm, &(req.back()));
	outgoing_ptr = outgoing_buf.size();

      }
      
      // Return 0 so that we don't print anything
      return 0;
    }
  }
}
#endif	

