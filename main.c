#include <stdio.h>
#include <stdlib.h>

#define LAB4_EXTEND
#include "Lab4_IO.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

void initializeRanks(double *r, int nodeCount);

int main() {
	int nodeCount, i, j;
	int *numInLinks, *numOutLinks;
	double *r, *rPre;
	double dampConst;

	struct node *nodeHead;

	// numInLinks: links that point to one node
	// numOutLinks: links that a node point to
    get_node_stat(&nodeCount, &numInLinks, &numOutLinks);

    printf("Node Count: %d\n", nodeCount);

    // Step 1: initialize all nodes to have same rank
    r = malloc(nodeCount*sizeof(double));
 	initializeRanks(r, nodeCount);

    // Step 2: Equation 3, 
    node_init(&nodeHead, numInLinks, numOutLinks, 0, nodeCount);
    
    // get our previous R array to know when to stop the process
    rPre = malloc(nodeCount*sizeof(double));
    
    int iterationcount = 0;
    dampConst = (1.0 - DAMPING_FACTOR) / nodeCount;
    do {
    	++iterationcount;
        vec_cp(r, rPre, nodeCount);
        for ( i = 0; i < nodeCount; ++i){
            r[i] = 0;
            for ( j = 0; j < nodeHead[i].num_in_links; ++j)
                r[i] += rPre[nodeHead[i].inlinks[j]] / numOutLinks[nodeHead[i].inlinks[j]];
            r[i] *= DAMPING_FACTOR;
            r[i] += dampConst;
        }
    } while(rel_error(r, rPre, nodeCount) >= EPSILON);

    Lab4_saveoutput(r, nodeCount, 100.0);
	return 0;
}

void initializeRanks(double *r, int nodeCount) {
   int i;
	for(i = 0; i < nodeCount; ++i)
		r[i] = (double)(1.0/nodeCount);
}