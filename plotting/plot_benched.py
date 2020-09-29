import os
import numpy as np
import matplotlib.pyplot as plt

plt_format = '.png'
bench_path = "../build/benchmark_data/"
alg_names  = ["grover"]

replot_all = False
plot_concur_perf_bool = True
plot_peaknodes_bool   = True
plot_histories_bool   = True


def plot_history(input_path, output_path, alg_name):
    ys = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)
    x  = np.arange(start=0, stop=len(ys), step=1)

    _, ax1 = plt.subplots()
    color1 = 'tab:orange'
    color2 = 'tab:blue'

    # gates vs nodes
    ax1.set_xlabel('gates')
    ax1.set_ylabel('qdd nodes', color=color1)
    ax1.plot(x, ys['nodes'], color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    # gates vs (total) ctable entries
    ax2 = ax1.twinx()
    ax2.set_ylabel('complex table entries', color=color2)
    ax2.plot(x, ys['amps'], color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()
    plt.close()


def plot_histories(histories_path, output_folder, alg_name):
    # make sure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # iterate over all .csv files
    for filename in os.listdir(histories_path):
        if filename.endswith(".csv"):
            hist_file_path = histories_path + filename
            output_path = output_folder + filename[:-4] + plt_format
            plot_history(hist_file_path, output_path, alg_name)


def plot_qubits_vs_peak_nodes(data_folder, alg_name):
    input_path  = data_folder + 'summary.csv'
    output_path = data_folder + 'peak_nodes' + plt_format
    data = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)
    x = data['qubits']
    y = data['peak_nodes']
    plt.scatter(x, y)
    plt.xlabel('qubits')
    plt.ylabel('qdd peak nodes')
    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()
    plt.close()

def plot_concurrency_performance(data_folder, alg_name):
    input_path  = data_folder + 'summary.csv'
    output_path = data_folder + 'concurrency' + plt_format
    data = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)

    lengend_entries = []
    
    # for different number of qubits
    qubits  = np.unique(data['qubits'])
    workers = np.unique(data['workers'])
    for q in qubits:
        subset = data[np.where(data['qubits'] == q)]
        speedups = np.zeros(workers.shape)
        
        # asume grover, group by flag
        flags = np.unique(subset['flag'])
        for flag in flags:
            subsubset = subset[np.where(subset['flag'] == flag)]
            speed_w1 = np.mean(subsubset[np.where(subsubset['workers'] == 1)]['runtime'])
            for i, w in enumerate(workers):
                speed_w = np.mean(subsubset[np.where(subsubset['workers'] == w)]['runtime'])
                speedups[i] += (speed_w / speed_w1)**(-1)

        # these are the speedups averaged over the different flags
        speedups /= flags.shape

        # actually plot stuff
        plt.scatter(workers, speedups)
        w1_times = np.round(subset[np.where(subset['workers'] == 1)]['runtime'], 3)
        leg = '{} qubits, time $w_1$ ({},{})'.format(int(q), np.min(w1_times), np.max(w1_times))
        lengend_entries.append(leg)

    plt.ylabel('average speedup')
    plt.xlabel('number of workers')
    plt.xticks(workers.astype(int))
    plt.legend(lengend_entries)
    plt.plot(workers, np.ones(workers.shape), color='grey', linestyle='--')
    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.savefig(output_path)
    plt.clf()
    plt.close()


# iterates over all folders in the bench_path and plots everything it can plot
def plot_all():

    # iterate over all algorithms
    for alg_name in alg_names:
        alg_path = bench_path + alg_name + "/"

        # iterate over all experiments
        for exp_folder in os.listdir(alg_path):
            exp_path = alg_path + exp_folder + "/"
            if ('concurrency.png' in os.listdir(exp_path)):
                print("skipping {}".format(exp_folder))
            else:
                print("plotting {}".format(exp_folder))

                # plot qubits vs peak nodes
                if (plot_peaknodes_bool):
                    plot_qubits_vs_peak_nodes(exp_path, alg_name)

                # plot concurrency performance
                if (plot_concur_perf_bool):
                    plot_concurrency_performance(exp_path, alg_name)

                # plot histories for all runs
                if (plot_histories_bool):
                    histories_path = exp_path + "run_histories/"
                    output_folder  = histories_path + "plots/"
                    plot_histories(histories_path, output_folder, alg_name)


if __name__ == "__main__":
    plot_all()
