import os
from datetime import datetime
import subprocess
import json
import re as regular_expression
import pandas as pd
import matplotlib.pyplot as plt

JSON_DIR = 'benchmark/json/'
#
# Function to flatten nested JSON
#
def flatten_json(nested_json, parent_key='', sep='.'):
    items = []
    for k, v in nested_json.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_json(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

#
# Process all json files into one dataframe.
#

def convert_json_file_to_json_data(
    json_filename : str, 
    timestamp_day : int, 
    timestamp_yymmdd : int, 
    softwareversion : str, 
    precision : int, 
    method : str):

    filepath_json = os.path.join(JSON_DIR, json_filename)
    print(filepath_json)

    # Convert json to dataframe
    with open(filepath_json, 'r') as file:
        data = json.load(file)

    data = flatten_json(data)

    # Extend the json data with file parameters    
    data["timestamp_day"] = timestamp_day
    data["timestamp_yymmdd"] = timestamp_yymmdd
    data["softwareversion"] = softwareversion
    data["precision"] = precision
    data["method"] = method

    return data


def convert_all_json_files_to_dataframe():

    # Capture all files in the directory
    all_files = os.listdir(JSON_DIR)

    # Select only the sql files
    json_files = [file for file in all_files if file.endswith('.json')]

    # Remove something if needed
    json_files = [file for file in json_files if not file.startswith('2022', 0)]

    # Sort the files
    json_files = sorted(json_files)

    # Add every json file into a dataframe row
    first_file = True

    for file in json_files:

        print(f"Processing json file: {file}")
        
        words = file.split('_')
        timestamp_day = words[0]
        timestamp_yymmdd = words[1]
        softwareversion = words[2]
        precision = words[3]
        method = words[4]

        json_data = convert_json_file_to_json_data(file, timestamp_day, timestamp_yymmdd, softwareversion, precision, method)

        # Append new row to dataframe
        if first_file:
            df = pd.DataFrame([json_data])
            first_file = False
        else:
            df.loc[len(df)] = json_data

    return df

#
# Plot functions
#

def generate_time_nr_qubits_plots(df):

    # Plotting time and nr qubits
    df_ae = (
        df[
            (df['statistics.benchmark'].str[0:2] == 'ae') & 
            (df['precision'] == '64') & 
            (df['method'] == 'MTBDD')
          ]
        .sort_values(by='statistics.n_qubits')
        .plot(
            x='statistics.n_qubits', 
            y='statistics.simulation_time', 
            kind='line', 
            marker='o', 
            color='blue')
    )
    df_ws = (
        df[
            (df['statistics.benchmark'].str[0:6] == 'wstate') & 
            (df['precision'] == '64') & 
            (df['method'] == 'MTBDD')
          ]
        .sort_values(by='statistics.n_qubits')
        .plot(
            x='statistics.n_qubits', 
            y='statistics.simulation_time', 
            kind='line', 
            marker='o', 
            color='red', 
            ax=plt.gca())
    )

    # 'r' (red)
    # 'g' (green)
    # 'b' (blue)
    # 'c' (cyan)
    # 'm' (magenta)
    # 'y' (yellow)
    # 'k' (black)
    # 'w' (white)

    # Display the plot
    plt.title('wstate red, ae blue')
    plt.xlabel('qubits')
    plt.ylabel('wall time (s)')
    plt.grid(True)
    plt.savefig('benchmark/plot/time_nr_qubits.pdf')

    return

#
# Time against time plots
#

def generate_time_time_method_plots(df):

    # Plotting time and nr qubits
    df_x = (
        df[
            (df['precision'] == '64') & 
            (df['method'] == 'QMDD')
          ]
    )
    df_y = (
        df[
            (df['precision'] == '64') & 
            (df['method'] == 'MTBDD')
          ]
    )

    # Merge the two DataFrames on the column to compare
    merged_df = df_x.merge(df_y, on='statistics.benchmark', how='outer', indicator=True)

    # Identify unequal values
    unequal_values = merged_df[merged_df['_merge'] != 'both']

    # Display the unequal values
    print("Unequal values:")
    print(unequal_values)

    fig, ax = plt.subplots()

    ax.scatter(
        df_x['statistics.simulation_time'], 
        df_y['statistics.simulation_time'], 
        color='blue',
        marker='o'
        )

    # Set logarithmic scale for both axes
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set limits for both axes
    ax.set_xlim(0.01, 1000)
    ax.set_ylim(0.01, 1000)

    # Display the plot
    plt.title('Wall time of all benchmarks, precision=64')
    plt.xlabel('QMDD wall time (s)')
    plt.ylabel('MTBDD wall time (s)')
    plt.grid(True)
    plt.savefig('benchmark/plot/time_time_method.pdf')

    return

#
# Time against time plots
#

def generate_time_time_per_node_plots(df):

    # Plotting time and nr qubits
    df_x = (
        df[
            (df['precision'] == '64') & 
            (df['method'] == 'QMDD')
          ]
    )
    df_x = df_x.assign(time_per_node = df_x['statistics.simulation_time'] / df_x['statistics.final_nodes'])

    df_y = (
        df[
            (df['precision'] == '64') & 
            (df['method'] == 'MTBDD')
          ]
    )
    df_y = df_y.assign(time_per_node = df_y['statistics.simulation_time'] / df_y['statistics.final_nodes'])

    fig, ax = plt.subplots()

    ax.scatter(
        df_x['time_per_node'], 
        df_y['time_per_node'], 
        color='blue',
        marker='o'
        )

    # Set logarithmic scale for both axes
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set limits for both axes
    ax.set_xlim(0.00001, 1)
    ax.set_ylim(0.00001, 1)

    # Display the plot
    plt.title('Wall time per node of all benchmarks, precision=64')
    plt.xlabel('QMDD wall time (s)')
    plt.ylabel('MTBDD wall time (s)')
    plt.grid(True)
    plt.savefig('benchmark/plot/time_time_per_node.pdf')

    return

#
# Time against time for different precisions plots
#

def generate_time_time_precision_plots(df):

    # Plotting time and nr qubits
    df_x = (
        df[
            (df['precision'] == '256') & 
            (df['method'] == 'MTBDD')
          ]
    )
    df_y = (
        df[
            (df['precision'] == '64') & 
            (df['method'] == 'MTBDD')
          ]
    )

    print(df_x)
    print(df_y)

    # Merge the two DataFrames on the column to compare
    merged_df = df_x.merge(df_y, on='statistics.benchmark', how='outer', indicator=True)

    # Identify unequal values
    unequal_values = merged_df[merged_df['_merge'] != 'both']

    # Display the unequal values
    print("Unequal values:")
    print(unequal_values)

    fig, ax = plt.subplots()

    ax.scatter(
        df_x['statistics.simulation_time'], 
        df_y['statistics.simulation_time'], 
        color='blue',
        marker='o'
        )

    # Set logarithmic scale for both axes
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Set limits for both axes
    ax.set_xlim(0.01, 1000)
    ax.set_ylim(0.01, 1000)

    # Display the plot
    plt.title('Wall time of all benchmarks, precision 64 against 256')
    plt.xlabel('MTBDD 256 wall time (s)')
    plt.ylabel('MTBDD 64 wall time (s)')
    plt.grid(True)
    plt.savefig('benchmark/plot/time_time_precision.pdf')

    return

#
# Main
#

df = convert_all_json_files_to_dataframe()

print(df)

generate_time_nr_qubits_plots(df)

generate_time_time_method_plots(df)

generate_time_time_per_node_plots(df)

generate_time_time_precision_plots(df)

#generate_minimum_precision_nr_qubits_plots()     # f(precision, nr_qubits) for one circuittype. method = MTBDD

#df.save()
