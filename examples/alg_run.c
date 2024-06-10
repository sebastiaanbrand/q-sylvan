#include <argp.h>
#include <inttypes.h>
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
static int rseed = 0;
static int depth = 0; // must be set for supremacy
static size_t min_tablesize = 1LL<<25;
static size_t max_tablesize = 1LL<<25;
static size_t min_cachesize = 1LL<<20;
static size_t max_cachesize = 1LL<<20;
static size_t wgt_tab_size  = 1LL<<23;
static double tolerance     = 1e-14;
static int wgt_table_type   = COMP_HASHMAP;
static int wgt_norm_strat   = NORM_MAX;
static int wgt_inv_caching  = 1;

static int grover_flag = 1; // 0 = random, 1 = 11..1

static int shor_N = 0;
static int shor_a = 0;

static char* csv_outputfile = NULL;

enum algorithms {
    alg_grover,
    alg_shor,
    alg_supremacy
};

static struct argp_option options[] =
{
    {"workers", 'w', "<workers>", 0, "Number of workers/threads (default=1)", 0},
    {"qubits", 'q', "<nqubits>", 0, "Number of qubits (must be set for Grover)", 0},
    {"rseed", 'r', "<random-seed>", 0, "Set random seed", 0},
    {"depth", 'd', "<depth>", 0, "Depth of circuits with arbitrary depth (e.g. supremacy)", 0},
    {"norm-strat", 's', "<low|max|min|l2>", 0, "Edge weight normalization strategy", 0},
    {"tol", 1, "<tolerance>", 0, "Tolerance for deciding edge weights equal (default=1e-14)", 0},
    {"inv-caching", 2, "<0|1>", 0, "Turn inverse chaching of edge weight computations on/off (default=on)", 0},
    {"grover-flag", 20, "<random|ones>", 0, "Grover flag (default=11..1)", 0},
    {"shor-N", 30, "<N>", 0, "N to factor with Shor's algorithm", 0},
    {"shor-a", 31, "<a>", 0, "value 'a' to use in Shor's algorithm (chosen random if not set)", 0},
    {"csv-output", 40, "<filename>", 0, "Write stats to given filename (or append if exists)", 0},
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
    case 'r':
        rseed = atoi(arg);
        break;
    case 'd':
        depth = atoi(arg);
        break;
    case 's':
        if (strcmp(arg, "low")==0) wgt_norm_strat = NORM_LOW;
        else if (strcmp(arg, "max")==0) wgt_norm_strat = NORM_MAX;
        else if (strcmp(arg, "min")==0) wgt_norm_strat = NORM_MIN;
        else if (strcasecmp(arg, "l2")==0) wgt_norm_strat = NORM_L2;
        else argp_usage(state);
        break;
    case 1:
        tolerance = atof(arg);
        break;
    case 2:
        wgt_inv_caching = atoi(arg);
        break;
    case 20:
        if (strcmp(arg, "random")==0) grover_flag = 0;
        else if (strcmp(arg, "ones")==0) grover_flag = 1;
        else argp_usage(state);
        break;
    case 30:
        shor_N = atoi(arg);
        break;
    case 31:
        shor_a = atoi(arg);
        break;
    case 40:
        csv_outputfile = arg;
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
    int success;
    int nqubits;
    QMDD final_qmdd;
    int64_t final_nodecount;
    double final_magnitude;
} stats_t;
stats_t stats = {0};

static void
write_csv_stats()
{
    FILE *fp = fopen(csv_outputfile, "a");
    // write header if file is empty
    fseek (fp, 0, SEEK_END);
        long size = ftell(fp);
        if (size == 0)
            fprintf(fp, "%s\n", "algorithm, nqubits, tolerance, norm-strat, inv-cache, workers, success, runtime, final_nodecount, final_magnitude");
    // append stats of this run
    int max_length = 100;
    char alg_name[max_length];
    if (algorithm == alg_grover) snprintf(alg_name, max_length, "grover");
    else if (algorithm == alg_shor) snprintf(alg_name, max_length, "shor");
    else if (algorithm == alg_supremacy) snprintf(alg_name, max_length, "supremacy-depth%d", depth);
    fprintf(fp, "%s, %d, %.3e, %d, %d, %d, %d, %lf, %" PRIu64 ", %0.5lf\n",
            alg_name,
            stats.nqubits,
            tolerance,
            wgt_norm_strat,
            wgt_inv_caching,
            workers,
            stats.success,
            stats.runtime,
            stats.final_nodecount,
            stats.final_magnitude);
    fclose(fp);
}

/*****************************</Info and logging>******************************/





/*******************************<Run algorithms>*******************************/

void
run_grover()
{
    if (qubits <= 0) Abort("--qubits=<num> must be set for Grover\n");
    stats.nqubits = qubits + 1;

    bool *flag;
    if (grover_flag == 1) flag = qmdd_grover_ones_flag(qubits+1);
    else flag = qmdd_grover_random_flag(qubits+1);

    // Run + time Grover
    INFO("Running Grover for %d qubits (+1 ancilla)\n", qubits);
    double t1 = wctime();
    stats.final_qmdd = qmdd_grover(qubits, flag);
    double t2 = wctime();
    stats.runtime = t2-t1;

    // Sanity checks on final state
    // 1. Check flag probability (marginalize ancilla qubit out)
    flag[qubits] = 0; AADD_WGT amp0 = aadd_getvalue(stats.final_qmdd, flag);
    flag[qubits] = 1; AADD_WGT amp1 = aadd_getvalue(stats.final_qmdd, flag);
    double flag_prob = qmdd_amp_to_prob(amp0) + qmdd_amp_to_prob(amp1);
    INFO("Measure Grover flag with prob %lf\n", flag_prob);

    if (flag_prob > 0.9 && flag_prob < 1.0) {
        stats.success = 1;
    }

    free(flag);
    INFO("Grover Time: %f\n", stats.runtime);
}

void
run_supremacy()
{
    if (depth <= 0) Abort("--depth=<depth> must be set for Supremacy\n");
    if (qubits != 5 && qubits != 20) Abort("--qubits=<5|20> for supremacy\n");
    stats.nqubits = qubits;

    double t1 = wctime();
    if (qubits == 5)        stats.final_qmdd = supremacy_5_1_circuit(depth);
    else if (qubits == 20)  stats.final_qmdd = supremacy_5_4_circuit(depth);
    double t2 = wctime();
    stats.runtime = t2-t1;

    // don't have a sanity check other than the magnitude of the final qmdd
    stats.success = -1;

    INFO("Supremacy-%d Time: %f\n", qubits, stats.runtime);
}

void
run_shor()
{
    if (shor_N <= 0) Abort("--shor-N=<N> must be set for Shor\n");
    stats.nqubits = shor_get_nqubits(shor_N);
    if (shor_a == 0) shor_a = shor_generate_a(shor_N);

    INFO("Running Shor with %d qubits to factor %d, with --shor-a=%d\n", stats.nqubits, shor_N, shor_a);
    double t1 = wctime();
    int factor = shor_run(shor_N, shor_a, false);
    double t2 = wctime();
    stats.runtime = t2-t1;
    stats.final_qmdd = shor_get_final_qmdd();

    if (factor != 0 && shor_N % factor == 0) {
        INFO("Shor: Found factor %d of %d\n", factor, shor_N);
        stats.success = 1;
    } else {
        INFO("Shor: Did not find factor\n");
    }
    
    INFO("Shor Time: %f\n", stats.runtime);
}

/******************************</Run algorithms>*******************************/





int main(int argc, char **argv)
{
    /* Parse arguments, set startup time for INFO messages. */
    argp_parse(&argp, argc, argv, 0, 0, 0);
    t_start = wctime();

    /* Init Lace + Sylvan */
    lace_start(workers, 0);
    sylvan_set_sizes(min_tablesize, max_tablesize, min_cachesize, max_cachesize);
    sylvan_init_package();
    qsylvan_init_simulator(wgt_tab_size, wgt_tab_size, tolerance, wgt_table_type, wgt_norm_strat);
    wgt_set_inverse_chaching(wgt_inv_caching);

    // e.g. for choosing random 'a' in Shor
    if (rseed == 0) rseed = time(NULL);
    srand(rseed);

    /* Print some info */
    INFO("Edge weight normalization: %d\n", wgt_norm_strat);
    INFO("Edge weight tolerance: %.3e\n", tolerance);
    INFO("Edge weight inverse caching: %d\n", wgt_inv_caching);
    INFO("Workers: %d\n", workers);
    INFO("Random seed: %d\n", rseed);

    /* Run the given quantum algorithm */
    if (algorithm == alg_grover) {
        run_grover();
    } else if (algorithm == alg_shor) {
        run_shor();
    } else if (algorithm == alg_supremacy) {
        run_supremacy();
    }

    /* Some stats */
    stats.final_magnitude = qmdd_get_magnitude(stats.final_qmdd, stats.nqubits);
    INFO("Magnitude of final state: %.05lf\n", stats.final_magnitude);
    stats.final_nodecount = aadd_countnodes(stats.final_qmdd);
    INFO("Final Nodecount: %" PRIu64 "\n", stats.final_nodecount);
    if (csv_outputfile != NULL) {
        INFO("Writing csv output to %s\n", csv_outputfile);
        write_csv_stats();
    }

    /* Cleanup */
    sylvan_quit();
    lace_stop();

    return 0;   
}
