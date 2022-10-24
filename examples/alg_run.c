#include <argp.h>
#include <stdlib.h>
#include <string.h>

#include "grover.h"
#include "shor.h"
#include "supremacy.h"


/**********************<Arguments (configured via argp)>***********************/

static int algorithm = 0;
static int qubits = 0; // must be set for Grover
static int workers = 1;
static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<20;
static size_t max_cachesize = 1LL<<20;
static size_t wgt_tab_size  = 1LL<<23;
static double tolerance     = 1e-14;
static int wgt_table_type   = COMP_HASHMAP;
static int wgt_norm_strat   = NORM_LARGEST;

enum algorithms {
    alg_grover,
    alg_shor,
    alg_supremacy
};

static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers (default=1)", 0},
    {"qubits", 'q', "<nqubits>", 0, "Number of qubits (must be set for Grover)", 0},
    {0, 0, 0, 0, 0, 0}
};
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
    switch (key) {
    case 'w':
        workers = atoi(arg);
        break;
    case 'q':
        qubits = atoi(arg);
        break;
    case ARGP_KEY_ARG:
        if (state->arg_num >= 1) argp_usage(state);
        else if (strcmp(arg, "grover")==0) algorithm = alg_grover;
        else if (strcmp(arg, "shor")==0) algorithm = alg_shor;
        else if (strcmp(arg, "supremacy")==0) {
            algorithm = alg_supremacy;
            if (qubits == 0) qubits = 20;
        }
        else argp_usage(state);
        break;
    case ARGP_KEY_END:
        if (state->arg_num < 1) argp_usage(state);
        break;
    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
static struct argp argp = { options, parse_opt, "<alg_name>", 0, 0, 0, 0 };

/*********************</Arguments (configured via argp)>***********************/


int main(int argc, char **argv)
{
    /* Parse arguments */
    argp_parse(&argp, argc, argv, 0, 0, 0);


    /* Init Lace + Sylvan */
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, tolerance, wgt_table_type, wgt_norm_strat);


    printf("TODO: run alg %d from command line\n", algorithm);

    /* Run the given quantum algorithm */
    if (algorithm == alg_grover) {
        // TODO: alg + timing
        printf("(WIP) run Grover on %d qubits\n", qubits);
    } else if (algorithm == alg_shor) {
        // TODO
        printf("(WIP) run Shor\n");
    } else if (algorithm == alg_supremacy) {
        // TODO
        printf("(WIP) run %d-qubit supremacy circuit\n", qubits);
    }


    /* Cleanup */
    sylvan_quit();
    lace_exit();

    return 0;   
}