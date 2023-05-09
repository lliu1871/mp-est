#include "neighbour_joining.h"
#include "tool.h"
#include "mpest.h"

/* Get the size in bytes of a dist_matrix with the given number of species */
size_t          dist_matrix_size (uint32_t species_count);
void            dist_matrix_free (dist_matrix *dmat);
double          dist_matrix_distance (const dist_matrix *dmat, uint32_t s1, uint32_t s2);
double          dist_matrix_avg_distance_from_others (const dist_matrix *dmat, uint32_t s);
void            dist_matrix_compute_avg_distances (const dist_matrix *dmat, double distances[]);
void            dist_matrix_print (const dist_matrix *dmat);
dist_matrix     *load_file(const char *file_name);
char            *neigh_strdup(const char *s);
size_t          trim_trailing_space(char *s);
size_t          filename_copy(const char *path, char *dest, size_t size);
size_t          btree_storage_size(uint32_t nodes_count);
btree_storage   *btree_storage_init(uint32_t nodes_count);
void            btree_storage_free(btree_storage *tree);
btree_node      *btree_storage_fetch(btree_storage *storage);
char            *btree_node_set_name(btree_node *node, const char *name);
uint32_t        btree_get_height(btree_node *root);
void            btree_print_tree(btree_node *root);
void            btree_print_trees(btree_node **trees, uint32_t tree_count);
void            nj_find_nearest_clusters(const dist_matrix *dmat, const double u[], uint32_t *c1, uint32_t *c2);
dist_matrix     *nj_join_clusters(const dist_matrix *dmat, const char *new_name, uint32_t c1, uint32_t c2);
btree_storage   *nj_tree_init(const dist_matrix *dmat, btree_node **leafs);
void            nj_tree_add_node(const dist_matrix *dmat, const double u[], btree_storage *storage, btree_node **partial_trees, const char *name, uint32_t c1, uint32_t c2);
//static void     btree_write_node_to_file (FILE *fout, btree_node *node);
static void     btree_write_node_to_string (btree_node *node);

//char njtreeOutfile[] = "njstree.tre";

int NeighborJoining (dist_matrix *dmat) {
        
        //size_t length = filename_copy(input_file, NULL, 0);
       // char filename[length + 1];
        //filename_copy(input_file, filename, sizeof(filename));
        
        //dist_matrix *dmat = load_file(input_file);
    
        //FILE *fout = fopen(njtreeOutfile, "w");
        
        if (!dmat) {
            // Error, stop here
            return EXIT_FAILURE;
        }
               
        double u[dmat->species_count];
        
        /* Ensure that cluster_name has enough space for the longest possible name */
        char cluster_name[2 + 3 * sizeof(dmat->species_count) + 1];
        uint32_t cluster_id = 1;
        
        btree_node *partial_trees[dmat->species_count];
        btree_storage *tree_storage = nj_tree_init(dmat, partial_trees);
        
        while (dmat->species_count >= 2) {
            
            /* Compute the average distance of each clusters from the others */
            dist_matrix_compute_avg_distances(dmat, u);
            
            /* Find the pair of nearest clusters */
            uint32_t c1, c2;
            nj_find_nearest_clusters(dmat, u, &c1, &c2);
            
            /* Generate a name for the new cluster */
            unsigned long result = snprintf(cluster_name, sizeof(cluster_name), "C_%" PRIu32, cluster_id);
            assert(result > 0 && result < sizeof(cluster_name));
            
            /* Add a node for the new cluster to the array of partial trees */
            nj_tree_add_node(dmat, u, tree_storage, partial_trees, cluster_name, c1, c2);
            
            /* Create a new dist_matrix joining the specified clusters */
            dist_matrix *joined = nj_join_clusters(dmat, cluster_name, c1, c2);
            
            if (joined == NULL) {
                /* Error, stop here */
                break;
            }
            
            /* Release the old distance matrix */
            dist_matrix_free(dmat);
            dmat = joined;
            
            cluster_id++;
        }
        
        btree_node *phyl_tree = partial_trees[0];
       
        /*print tree on screen*/
        //btree_print_tree(phyl_tree);
        //printf("\n\n");

        /*save tree to file*/
        //btree_write_node_to_file(fout, phyl_tree);
        //fprintf(fout,";");
        //fclose(fout);

        /*save tree to string*/
        printStringSize = 200;
        printString = (char *)malloc(printStringSize * sizeof(char));
        if (!printString){
            Print ("Problem allocating printString (%d)\n", printStringSize * sizeof(char));
            return ERROR;
        }
        *printString = '\0';
        btree_write_node_to_string (phyl_tree);
        strcat(printString,";");
        
        dist_matrix_free(dmat);
        btree_storage_free(tree_storage);
             
    return EXIT_SUCCESS;
}

/*static void btree_write_node_to_file (FILE *fout, btree_node *node) {
    if (node->left->left == NULL) {
        fprintf(fout, "(%s:%lf,", node->left->node_name, node->distance_left);
    }else{
        fprintf(fout, "(");
        btree_write_node_to_file(fout, node->left);
        fprintf(fout, ":%lf,", node->distance_left);
    }

    if (node->right->left == NULL) {
        fprintf(fout, "%s:%lf)", node->right->node_name, node->distance_right);
    }else{
        btree_write_node_to_file(fout, node->right);
        fprintf(fout, ":%lf)", node->distance_right);
    }
}*/

static void btree_write_node_to_string(btree_node *node){
    char	*tempStr;
	int     tempStrSize = 200;

	tempStr = (char *) malloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr){
		Print ("Problem allocating tempString (%d)\n", tempStrSize * sizeof(char));
	}
	
    if (node->left->left == NULL) {
        SaveSprintf (&tempStr, &tempStrSize, "(%s:%lf,", node->left->node_name, node->distance_left);
        AddToPrintString (tempStr);
    }else{
        SaveSprintf (&tempStr, &tempStrSize, "(");
		AddToPrintString (tempStr);
        btree_write_node_to_string(node->left);
        SaveSprintf (&tempStr, &tempStrSize, ":%lf,", node->distance_left);
        AddToPrintString (tempStr);
    }

    if (node->right->left == NULL) {
        SaveSprintf (&tempStr, &tempStrSize, "%s:%lf)", node->right->node_name, node->distance_right);
        AddToPrintString (tempStr);
    }else{
        btree_write_node_to_string(node->right);
        SaveSprintf (&tempStr, &tempStrSize, ":%lf)", node->distance_right);
        AddToPrintString (tempStr);
    }
}

void nj_find_nearest_clusters(const dist_matrix *dmat, const double u[], uint32_t *c1, uint32_t *c2) {
    assert(dmat->species_count >= 2);
    assert(c1 != NULL);
    assert(c2 != NULL);
   
    /*
     * Find the two clusters that minimize the distance metric:
     *      min { D(c1, c2) - U(c1) - U(c2) } among all c1, c2 in dmat
     */
    double min_distance = INFINITY;

    for (uint32_t i = 0; i < dmat->species_count; i++) {
        for (uint32_t j = 0; j < i; j++) {
            double distance = dist_matrix_distance(dmat, i, j) - u[i] - u[j];

            if (distance < min_distance) {
                min_distance = distance;

                *c1 = i;
                *c2 = j;
            }
        }
    }

    /* A pair of clusters should always be found */
    assert(isfinite(min_distance));
}

dist_matrix *nj_join_clusters(const dist_matrix *dmat, const char *new_name, uint32_t c1, uint32_t c2) {
    assert(dmat->species_count >= 2);

    dist_matrix *out = dist_matrix_init(dmat->species_count - 1);

    if (out == NULL) {
        perror("Unable to create output distance matrix");
        return NULL;
    }
    
    uint32_t k = 0;

    for (uint32_t i = 0; i < dmat->species_count; i++) {
        const char *species_name;
        uint32_t cluster_size;

        if (i == c1) {
            /* Replace c1 with the union of clusters c1 and c2 */
            species_name = new_name;
            cluster_size = dmat->cluster_sizes[c1] + dmat->cluster_sizes[c2];
        } else if (i == c2) {
            /* Remove cluster c2 */
            continue;
        } else {
            /* Leave other clusters unchanged */
            species_name = dmat->species_names[i];
            cluster_size = dmat->cluster_sizes[i];
        }

        dist_matrix_set_species_name(out, k, species_name);
        out->cluster_sizes[k] = cluster_size;

        /* Compute the distances */
        uint32_t l = 0;

        for (uint32_t j = 0; j < i; j++) {
            if (j == c2) {
                /* Remove cluster c2 */
                continue;
            }

            double distance;

            if (i == c1) {
                double d1j = dist_matrix_distance(dmat, c1, j);
                double d2j = dist_matrix_distance(dmat, c2, j);
                double d12 = dist_matrix_distance(dmat, c1, c2);

                distance = (d1j + d2j - d12) / 2;
            } else if (j == c1) {
                double di1 = dist_matrix_distance(dmat, i, c1);
                double di2 = dist_matrix_distance(dmat, i, c2);
                double d12 = dist_matrix_distance(dmat, c1, c2);

                distance = (di1 + di2 - d12) / 2;
            } else {
                distance = dist_matrix_distance(dmat, i, j);
            }

            *(dist_matrix_element(out, k, l)) = distance;

            l++;
        }

        k++;
    }

    return out;
}

btree_storage *nj_tree_init(const dist_matrix *dmat, btree_node **leafs) {
    uint32_t node_count = 2 * dmat->species_count - 1;
    btree_storage *storage = btree_storage_init(node_count);
    
    for (uint32_t i = 0; i < dmat->species_count; i++) {
        leafs[i] = btree_storage_fetch(storage);

        btree_node_set_name(leafs[i], dmat->species_names[i]);
    }
    
    return storage;
}

void nj_tree_add_node(const dist_matrix *dmat, const double u[], btree_storage *storage, btree_node **partial_trees, const char *name, uint32_t c1, uint32_t c2) {
    btree_node *node = btree_storage_fetch(storage);

    btree_node_set_name(node, name);

    node->left = partial_trees[c1];
    node->right = partial_trees[c2];
    
    double distance = dist_matrix_distance(dmat, c1, c2);
    node->distance_left = (distance + u[c1] - u[c2]) / 2;
    node->distance_right = (distance + u[c2] - u[c1]) / 2;

    partial_trees[c1] = node;

    for (uint32_t i = c2 + 1; i < dmat->species_count; i++) {
        partial_trees[i - 1] = partial_trees[i];
    }
}

/*************************************************************
    phylogenetic tree
**************************************************************/
size_t btree_storage_size(uint32_t nodes_count) {
    return sizeof(btree_storage) + (nodes_count * member_size(btree_storage, nodes[0]));
}

btree_storage *btree_storage_init(uint32_t nodes_count) {
    size_t size = btree_storage_size(nodes_count);

    /* Zero-initialize the storage to ensure that btree_storage_free can correctly
     * free the strings */
    btree_storage *storage = calloc(size, 1);

    if (storage != NULL) {
        storage->nodes_count = nodes_count;
    }
    
    return storage;
}

void btree_storage_free(btree_storage *storage) {
    if (storage) {
        for (uint32_t i = 0; i < storage->nodes_count; ++i) {
            free(storage->nodes[i].node_name);
        }

        free(storage);
    }
}

btree_node *btree_storage_fetch(btree_storage *storage) {
    assert(storage->used_nodes < storage->nodes_count);
    
    uint32_t idx = storage->used_nodes;
    btree_node *node = &storage->nodes[idx];
    
    storage->used_nodes++;
    
    return node;
}

char *btree_node_set_name(btree_node *node, const char *name) {
    free(node->node_name);

    node->node_name = neigh_strdup(name);

    return node->node_name;
}

uint32_t btree_get_height(btree_node *root) {
    if (root == NULL) {
        return 0;
    }
    
    uint32_t left = btree_get_height(root->left);
    uint32_t right = btree_get_height(root->right);

    uint32_t height = 1;
    height += (left >= right) ? left : right;
    
    return height;
}

static void btree_print_node(btree_node *root, double distance, uint32_t depth, bool *is_open) {
    if (root == NULL) {
        return;
    }

    btree_print_node(root->right, root->distance_right, depth + 1, is_open);

    if (depth > 0) {
        for (uint32_t i = 1; i < depth; i++) {
            printf("  %c        ", is_open[i] ? '|' : ' ');
        }
    
        printf("  |--%.2lf-- ", distance);
    }

    is_open[depth] = !is_open[depth];

    printf("%s\n", root->node_name);

    btree_print_node(root->left, root->distance_left, depth + 1, is_open);
}

void btree_print_tree(btree_node *root) {
    assert(root != NULL);
    
    uint32_t height = btree_get_height(root);
    
    bool is_open[height];
    
    for (uint32_t i = 0; i < height; i++) {
        is_open[i] = false;
    }
    
    btree_print_node(root, 0, 0, is_open);
}

void btree_print_trees(btree_node **trees, uint32_t tree_count) {
    for (uint32_t i = 0; i < tree_count; i++) {
        btree_print_tree(trees[i]);
        printf("\n");
    }
}

/*************************************************
           Utilities
*************************************************/
char *neigh_strdup(const char *src) {
    char *dst = NULL;

    if (src != NULL) {
        size_t length = strlen(src);

        dst = malloc(length + 1);

        if (dst != NULL) {
            strcpy(dst, src);
        }
    }

    return dst;
}

size_t trim_trailing_space(char *s) {
    size_t length = strlen(s);
    
    while (length > 0 && isspace(s[length - 1])) {
        length--;
    }
    
    s[length] = '\0';
    
    return length;
}

size_t filename_copy(const char *path, char *dest, size_t size) {
    size_t length = strlen(path);
    size_t end = length;
    
    while (end > 0 && path[end - 1] != '.') {
        end--;
    }
    
    if (end == 0) {
        /* No extension found */
        end = length;
    } else {
        /* Remove the trailing dot */
        end--;
    }
    
    size_t start = end;
    
    while (start > 0 && path[start - 1] != '/') {
        start--;
    }
    
    size_t count = end - start;
    
    if (dest != NULL) {
        size_t n = ((size - 1) < count) ? (size - 1) : count;
        
        memcpy(dest, path + start, n);
        dest[n] = '\0';
    }
    
    return count;
}

/***************************************************
    io
****************************************************/
#define CHECK_SCANF_RESULT(result, value, message, file, dmat) \
    if (result != value) {                                     \
        perror(message);                                       \
                                                               \
        dist_matrix_free(dmat);                                \
        fclose(file);                                          \
                                                               \
        return NULL;                                           \
    }

dist_matrix *load_file(const char *file_name) {
    FILE *f = fopen(file_name, "r");

    if (!f) {
        perror("Unable to open input file");
        return NULL;
    }

    int result;
    uint32_t species_count;

    result = fscanf(f, "%" SCNu32, &species_count);

    CHECK_SCANF_RESULT(result, 1, "Invalid species count", f, NULL);

    dist_matrix *dmat = dist_matrix_init(species_count);

    if (!dmat) {
        perror("Unable to create distance matrix");
        return NULL;
    }

    for (uint32_t i = 0; i < species_count; i++) {
        /* species name: up to 30 alphabetic or whitespace characters */
        char species_name[31];

        result = fscanf(f, " %30[^0-9*\n]", species_name);

        CHECK_SCANF_RESULT(result, 1, "Invalid species name", f, dmat);

        trim_trailing_space(species_name);

        dist_matrix_set_species_name(dmat, i, species_name);
        dmat->cluster_sizes[i] = 1;

        for (uint32_t j = 0; j < i; j++) {
            double *element = dist_matrix_element(dmat, i, j);

            result = fscanf(f, "%lf", element);

            CHECK_SCANF_RESULT(result, 1, "Invalid distance", f, dmat);
        }

        /* Ignore asterisks at the end of the line */
        //fscanf(f, " *");
    }

    fclose(f);

    return dmat;
}

/*****************************************
    distance
******************************************/

size_t dist_matrix_size(uint32_t species_count) {
    size_t matrix_size = species_count * (species_count - 1) / 2;
    size_t struct_size = sizeof(dist_matrix) + matrix_size * member_size(dist_matrix, distances[0]);

    return struct_size;
}

dist_matrix *dist_matrix_init(uint32_t species_count) {
    size_t size = dist_matrix_size(species_count);
    dist_matrix *dmat = malloc(size);

    if (dmat != NULL) {
        dmat->species_count = species_count;

        /* Zero-initialize the array of species names to ensure that dist_matrix_free
         * can correctly free the strings */
        dmat->species_names = calloc(species_count, member_size(dist_matrix, species_names[0]));

        if (dmat->species_names == NULL) {
            dist_matrix_free(dmat);

            return NULL;
        }
        
        dmat->cluster_sizes = malloc(species_count * member_size(dist_matrix, cluster_sizes[0]));
        
        if (dmat->cluster_sizes == NULL) {
            dist_matrix_free(dmat);
            
            return NULL;
        }
    }

    return dmat;
}

void dist_matrix_free(dist_matrix *dmat) {
    if (dmat) {
        if (dmat->species_names) {
            for (uint32_t i = 0; i < dmat->species_count; i++) {
                free(dmat->species_names[i]);
            }

            free(dmat->species_names);
        }
        
        if (dmat->cluster_sizes) {
            free(dmat->cluster_sizes);
        }

        free(dmat);
    }
}

char *dist_matrix_set_species_name(dist_matrix *dmat, uint32_t s, const char *species_name) {
    free(dmat->species_names[s]);

    dmat->species_names[s] = neigh_strdup(species_name);

    return dmat->species_names[s];
}

static uint32_t dist_matrix_get_offset(const dist_matrix *dmat, uint32_t s1, uint32_t s2) {
    assert(s1 != s2);
    assert(s1 < dmat->species_count);
    assert(s2 < dmat->species_count);

    if (s1 < s2) {
        return dist_matrix_get_offset(dmat, s2, s1);
    } else {
        return (s1 * (s1 - 1) / 2) + s2;
    }
}

double *dist_matrix_element(dist_matrix *dmat, uint32_t s1, uint32_t s2) {
    return dmat->distances + dist_matrix_get_offset(dmat, s1, s2);
}

double dist_matrix_distance(const dist_matrix *dmat, uint32_t s1, uint32_t s2) {
    return *(dmat->distances + dist_matrix_get_offset(dmat, s1, s2));
}

double dist_matrix_avg_distance_from_others(const dist_matrix *dmat, uint32_t s) {
    assert(dmat->species_count >= 2);
    
    double distance = 0;
    
    if (dmat->species_count > 2) {
        for (uint32_t i = 0; i < dmat->species_count; i++) {
            if (i == s) {
                continue;
            }
            
            distance += dist_matrix_distance(dmat, s, i);
        }

        distance /= (dmat->species_count - 2);
    }
    
    return distance;
}

void dist_matrix_compute_avg_distances(const dist_matrix *dmat, double distances[]) {
    for (uint32_t i = 0; i < dmat->species_count; i++) {
        distances[i] = dist_matrix_avg_distance_from_others(dmat, i);
    }
}

void dist_matrix_print(const dist_matrix *dmat) {
    
    /* Max value among species_name, the string "c_size" and the distances */
    int max_length = 6;

    for (uint32_t i = 0; i < dmat->species_count; i++) {
        int length = (int)strlen(dmat->species_names[i]);

        if (length > max_length) {
            max_length = length;
        }
    }

    printf("%*s\t", max_length, "");

    for (uint32_t i = 0; i < dmat->species_count; i++) {
        printf("%*s\t", max_length, dmat->species_names[i]);
    }

    printf("\n");

    for (uint32_t i = 0; i < dmat->species_count; i++) {
        printf("%*s\t", max_length, dmat->species_names[i]);

        for (uint32_t j = 0; j < i; j++) {
            printf("%*.4lf\t", max_length, dist_matrix_distance(dmat, i, j));
        }

        printf("%*s\n", max_length, "*");
    }

    printf("\n");
    printf("%*s\t", max_length, "c_size");
    
    
    for (uint32_t i = 0; i < dmat->species_count; i++) {
        printf("%*" PRIu32 "\t", max_length, dmat->cluster_sizes[i]);
    }
    
    printf("\n");

}
