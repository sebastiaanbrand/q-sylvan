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

#
# Main
#

df = convert_all_json_files_to_dataframe()

print(df)

#generate_time_nr_qubits_plots()                 # f(time(method), nr_qubit(circuittype)) for one method

#generate_time_time_plots()                      # f(time(EVDD), time(MTBDD)) for all circuits selected on softwareversion

#generate_time_precision_plots()                 # f(time(MTBDD), precision) for one circuittype, one nr of qubits

#generate_minimum_precision_nr_qubits_plots()    # f(precision, nr_qubits) for one circuittype. method = MTBDD

#df.save()
