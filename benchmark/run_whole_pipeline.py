#
# Run this pipeline in the <root>/benchmark directory!
#

import globals_pipeline as g

import generate_qasm_random_clif_T_gate as M_qasm
import convert_qasm_to_json as M_json
import convert_json_to_plots as M_plots

#
# QASM generation
#

# Generate and save circuits
for num_qubits in range(g.MIN_NUM_QUBITS, g.MAX_NUM_QUBITS, g.STEPSIZE_NUM_QUBITS):
    for num_gates in range(g.MIN_NUM_GATES, g.MAX_NUM_GATES, g.STEPSIZE_NUM_GATES):
        M_qasm.generate_and_save_circuits(num_qubits, num_gates)

#
# JSON generation
#
M_json.convert_all_qasm_files('v1.0.0')

#
# PLOTS generation
#
df = M_plots.convert_all_json_files_to_dataframe()

print("Conversion to flatten dataframe done")
#print(df)
#print(df.columns)

#print(df['statistics.norm'][df['statistics.norm'] < 1.0])

#generate_time_nr_qubits_plots(df)

#generate_time_time_method_plots(df)

#generate_time_time_per_node_plots(df)

#generate_time_time_precision_plots(df)

#generate_norm_norm_method_plots(df)

#generate_norm_norm_precision_plots(df)

M_plots.generate_time_nr_qubits_plots(df)

#df.save()





