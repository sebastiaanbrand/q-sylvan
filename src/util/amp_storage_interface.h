#ifndef AMP_STORAGE_INTERFACE
#define AMP_STORAGE_INTERFACE

typedef enum amp_storage_backend {
    COMP_HASHMAP, 
    REAL_HASHMAP, 
    REAL_TREE
} amp_storage_backend_t;

#define n_backends 3

// create(uint64_t size, double tolerance)
void * (*amp_store_create[n_backends])();

// free
void (*amp_store_free[n_backends])();

// get tolerance
double (*amp_store_get_tol[n_backends])();

void init_amp_storage_functions();

#endif // AMP_STORAGE_INTERFACE