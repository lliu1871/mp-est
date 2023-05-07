#define NEIGH_NEIGHBOUR_JOINING_H

#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <memory.h>
#include <inttypes.h>

#define member_size(type, member) sizeof(((type *)0)->member)

typedef struct dist_matrix dist_matrix;

struct dist_matrix {
    /* Number of species */
    uint32_t species_count;
    
    /* Names of the species */
    char **species_names;
    
    /* Number of species in each cluster */
    uint32_t *cluster_sizes;
    
    /* Distance matrix */
    double distances[];
};

typedef struct btree_node btree_node;

struct btree_node {
    /* Name of the node */
    char *node_name;

    /* Left and right children */
    btree_node *left, *right;
    
    /* Distances to the left and right children */
    double distance_left, distance_right;
};

typedef struct btree_storage btree_storage;

struct btree_storage {
    /* Number of nodes */
    uint32_t nodes_count;
    
    /* Number of used nodes / Index of first available node */
    uint32_t used_nodes;
    
    /* Array of tree nodes */
    btree_node nodes[];
};

/* wrapup function NJ*/
int NeighborJoining (dist_matrix *dmat);
dist_matrix *dist_matrix_init(uint32_t species_count);
char *dist_matrix_set_species_name(dist_matrix *dmat, uint32_t s, const char *species_name);
double *dist_matrix_element(dist_matrix *dmat, uint32_t s1, uint32_t s2);





