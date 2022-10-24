#include <argp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

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

static int grover_flag = 1; // 0 = random, 1 = 11..1

enum algorithms {
    alg_grover,
    alg_shor,
    alg_supremacy
};

static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers (default=1)", 0},
    {"qubits", 'q', "<nqubits>", 0, "Number of qubits (must be set for Grover)", 0},
    {"grover-flag", 20, "<random|ones>", 0, "Grover flag (default=11..1)", 0},
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
    case 20:
        if (strcmp(arg, "random")==0) grover_flag = 0;
        else if (strcmp(arg, "ones")==0) grover_flag = 1;
        else argp_usage(state);
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





/******************************<Info and logging>******************************/

/**
 * Obtain current wallclock time
 */
static double
wctime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

static double t_start;
#define INFO(s, ...) fprintf(stdout, "[% 8.2f] " s, wctime()-t_start, ##__VA_ARGS__)
#define Abort(...) { fprintf(stderr, __VA_ARGS__); fprintf(stderr, "Abort at line %d!\n", __LINE__); exit(-1); }

typedef struct stats {
    double runtime;
} stats_t;
stats_t stats = {0};

/*****************************</Info and logging>******************************/





/*******************************<Run algorithms>*******************************/

void
run_grover()
{
    if (qubits <= 0) {
        Abort("--qubits=<num> must be set for Grover\n");
    }

    bool *flag;
    if (grover_flag == 1) flag = qmdd_grover_ones_flag(qubits);
    else flag = qmdd_grover_random_flag(qubits);

    INFO("Running Grover for %d qubits\n", qubits);
    qmdd_grover(qubits, flag);

    free(flag);
}

/******************************</Run algorithms>*******************************/





int main(int argc, char **argv)
{
    /* Parse arguments, set startup time for INFO messages. */
    argp_parse(&argp, argc, argv, 0, 0, 0);
    t_start = wctime();


    /* Init Lace + Sylvan */
    lace_init(workers, 0);
    lace_startup(0, NULL, NULL);
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, tolerance, wgt_table_type, wgt_norm_strat);


    /* Run the given quantum algorithm */
    if (algorithm == alg_grover) {
        double t1 = wctime();
        run_grover();
        double t2 = wctime();
        stats.runtime = t2-t1;
        INFO("Grover Time: %f\n", stats.runtime);
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