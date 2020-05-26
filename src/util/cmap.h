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

static const double TOLERANCE = 1e10f; // (inverse of) float eq. tolerance

typedef struct complex_s {
  double r,i;
} complex_t;

/**
\typedef Lockless hastable database.
*/
typedef struct cmap_s cmap_t;

typedef size_t ref_t;

/**
\brief Create a new database.
\param len The length of the vectors to be stored here
\return the hashtable
*/
extern cmap_t *cmap_create (int size);

/**
\brief Free the memory used by a dbs.
*/
extern void cmap_free (cmap_t *dbs);

/**
\brief Find a vector with respect to a database and insert it if it cannot be fo
und.
\param dbs The dbs
\param vector The int vector
\retval idx The index that the vector was found or inserted at
\return 1 if the vector was present, 0 if it was added
*/
extern bool cmap_find_or_put (const cmap_t *dbs, const complex_t *v, ref_t *ret);

extern complex_t *cmap_get (const cmap_t *dbs, const ref_t ref);


#endif // CMAP
