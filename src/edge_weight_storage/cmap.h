#ifndef CMAP_H
#define CMAP_H

/**
\file cmap.h
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
#include "flt.h"


/**
\brief Create a new database.
\param len The length of the vectors to be stored here
\return the hashtable
*/
extern void *cmap_create(uint64_t size, double tolerance);

extern double cmap_get_tolerance();

/**
\brief Free the memory used by a dbs.
*/
extern void cmap_free(void *dbs);

/**
\brief Find a vector with respect to a database and insert it if it cannot be fo
und.
\param dbs The dbs
\param vector The int vector
\retval idx The index that the vector was found or inserted at
\return 1 if the vector was present, 0 if it was added, -1 if table was full
*/
extern int cmap_find_or_put(const void *dbs, const void *v, uint64_t *ret);

extern void * cmap_get(const void *dbs, const uint64_t ref);

extern uint64_t cmap_count_entries(const void *dbs);

extern void print_bitvalues(const void *dbs, const uint64_t ref);

#endif // CMAP
