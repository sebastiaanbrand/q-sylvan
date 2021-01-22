#ifndef RMAP_H
#define RMAP_H

/**
\file rmap.h
\brief Lockless non-resizing hash table implementation for fixed-length keys

@inproceedings{Laarman:2010:BMR:1998496.1998541,
  author = {Laarman, Alfons and van de Pol, Jaco and Weber, Michael},
  title = {{Boosting Multi-Core Reachability Performance with Shared Hash Tables}},
  booktitle = {Proceedings of the 2010 Conference on Formal Methods in Computer-Aided Design},
  series = {FMCAD '10},
  year = {2010},
  location = {Lugano, Switzerland},
  pages = {247--256},
  numpages = {10},
  url = {http://eprints.eemcs.utwente.nl/19281/},
  acmid = {1998541},
  publisher = {FMCAD Inc},
  address = {Austin, TX},
}
*/

#include <stdbool.h>
#include <stdint.h>
#include <string.h>


/**
\typedef Lockless hastable database.
*/
typedef struct rmap_s rmap_t;

typedef size_t ref_t;

/**
\brief Create a new database.
\param len The length of the vectors to be stored here
\return the hashtable
*/
extern rmap_t *rmap_create (uint64_t size, double tolerance);

extern double rmap_get_tolerance();

/**
\brief Free the memory used by a dbs.
*/
extern void rmap_free (rmap_t *dbs);

/**
\brief Find a vector with respect to a database and insert it if it cannot be fo
und.
\param dbs The dbs
\param vector The int vector
\retval idx The index that the vector was found or inserted at
\return 1 if the vector was present, 0 if it was added, -1 if table was full
*/
extern int rmap_find_or_put (const rmap_t *dbs, const double *v, ref_t *ret);

extern double *rmap_get (const rmap_t *dbs, const ref_t ref);

extern uint64_t rmap_count_entries (const rmap_t *rmap);

extern void rmap_print_bitvalues(const rmap_t *dbs, const ref_t ref);

#endif // RMAP
