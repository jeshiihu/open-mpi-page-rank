#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "timer.h"

#define LAB4_EXTEND
#include "Lab4_IO.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

void initializeRanks(double *r, int nodeCount);

int main() {
	int nodeCount, i, j;
	int *numInLinks, *numOutLinks;
	double *r, *rPre, *local_r;
	double dampConst;
	double begin, end, elapsed;
	int npes, myrank;
	MPI_Comm comm;

	struct node *nodeHead;

	// numInLinks: links that point to one node
	// numOutLinks: links that a node point to

	MPI_Init(NULL, NULL);
	comm = MPI_COMM_WORLD;
	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &myrank);

	get_node_stat(&nodeCount, &numInLinks, &numOutLinks);
	int lower = nodeCount*myrank / npes;
	int higher = nodeCount*(myrank+1) / npes;
	int local_node = nodeCount / npes;	

    // Step 1: initialize all nodes to have same rank
    r = malloc(nodeCount*sizeof(double));
 	initializeRanks(r, nodeCount);

	// Allocate memory for local r
	local_r = malloc(local_node*sizeof(double));

    // Step 2: Equation 3, 
    node_init(&nodeHead, numInLinks, numOutLinks, lower, higher);
    
    // get our previous R array to know when to stop the process
    rPre = malloc(nodeCount*sizeof(double));
    
	GET_TIME(begin);
    int iterationcount = 0;
    dampConst = (1.0 - DAMPING_FACTOR) / nodeCount;
    do {
	// Make a copy of our vector r
    	++iterationcount;
        vec_cp(r, rPre, nodeCount);


        for ( i = 0; i < local_node; ++i){
            local_r[i] = 0;
            for ( j = 0; j < nodeHead[i].num_in_links; ++j){

                local_r[i] += rPre[nodeHead[i].inlinks[j]] / numOutLinks[nodeHead[i].inlinks[j]];
			}
            local_r[i] *= DAMPING_FACTOR;
            local_r[i] += dampConst;
        }

		MPI_Allgather(local_r, local_node, MPI_DOUBLE, r, local_node, MPI_DOUBLE, comm);
    } while(rel_error(r, rPre, nodeCount) >= EPSILON);


	GET_TIME(end);
	elapsed = end - begin;
	MPI_Finalize();
	if (myrank == 0){
    		Lab4_saveoutput(r, nodeCount, elapsed);
	}

	node_destroy(nodeHead, local_node);
	free(r);
	free(rPre);
	free(local_r);
	return 0;
}

void initializeRanks(double *r, int nodeCount) {
   int i;
	for(i = 0; i < nodeCount; ++i)
		r[i] = (double)(1.0/nodeCount);
}
