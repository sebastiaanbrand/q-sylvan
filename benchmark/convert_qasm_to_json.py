import os
from datetime import datetime
import subprocess
import json
import re as regular_expression

SIM_QASM_EXE = './build/qasm/run_qasm_on_mtbdd'
#SIM_QASM_EXE = './build/qasm/run_qasm_on_qmdd'

QASM_DIR = 'benchmark/qasm/'
JSON_DIR = 'benchmark/json/'

#
# Simulate given quantum circuit in separate thread and store output in json.
#

def convert_qasm_to_json(qasm_filename : str, softwareversion : str):

    filepath_qasm = os.path.join(QASM_DIR, qasm_filename)
    print(filepath_qasm)

    # Split the filename into name and extension
    benchmarkname, ext = os.path.splitext(qasm_filename)
    if benchmarkname.split('_')[0] == 'grover-noancilla':  # out of memory for number 9
        return
    if benchmarkname.split('_')[0] == 'grover-v-chain':  # rccx currently unsupported
        return

    # Start simulating in separate thread
    output = subprocess.run(
        [SIM_QASM_EXE, filepath_qasm], #, '--state-vector'],
        stdout=subprocess.PIPE, check=False)

    # Convert output simulator to json format
    data = json.loads(output.stdout)

    # Get the current timestamp in the desired format (e.g., YYYYMMDD_HHMMSS)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create the new filename with the timestamp appended before the extension
    precision = 16
    method = "MTBDD"
    ext = ".json"
    json_filename = f"{timestamp}_{softwareversion}_{precision}_{method}_{benchmarkname}{ext}"

    filepath_json = os.path.join(JSON_DIR, json_filename)

    # Store json in file with timestamp
    with open(filepath_json, "w") as json_file:
        json.dump(data, json_file, indent=4)

    return


def convert_all_qasm_files(softwareversion : str):

    # Capture all files in the directory
    all_files = os.listdir(QASM_DIR)

    # Select only the sql files
    qasm_files = [file for file in all_files if file.endswith('.qasm')]

    # Remove something if needed
    qasm_files = [file for file in qasm_files if not file.startswith('2022', 0)]
    
    # Remove with regular expression, in this case: 
    #pattern = r"^(a|b|c|d|e|f|[g]).*" # remove files starting with a or b or ...
    #qasm_files = [file for file in qasm_files if not regular_expression.match(pattern, file)]

    #   runs into benchmark/qasm/grover-noancilla_indep_qiskit_10.qasm error
    #   GNU MP: Cannot allocate memory (size=17.565.568.910.312)

    #   runs into benchmark/qasm/*qaoa*.qasm error
    #   rzz unsupported gate 

    #   Processing qasm file takes too long: pricingcall_indep_qiskit_19.qasm
    #   benchmark/qasm/pricingcall_indep_qiskit_19.qasm

    # Sort the files
    qasm_files = sorted(qasm_files)

    for file in qasm_files:
        print(f"Processing qasm file: {file}")
        convert_qasm_to_json(file, softwareversion)
    return


#
# Main
#
convert_all_qasm_files('v1.0.0')

