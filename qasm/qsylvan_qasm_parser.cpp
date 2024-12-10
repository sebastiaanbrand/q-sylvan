/**
 * Copyright 2024 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#include "qsylvan_qasm_parser.h"
#include "parse_math/eval_expr.hpp"

/**
 * Minor shortcomings of this parser:
 * 
 * 1) Assumes only one instruction per line
 * 2) Ignores includes
 * 3) Assumes that for given rotation angles (e.g. rz(pi/4)), the expression 
 *    between the brackets doesn't contain additional brackets.
 * 
 */

std::vector<std::string> split(std::string to_split, std::string delims)
{
    std::vector<std::string> result;
    std::size_t index;

    while(to_split.size() > 0) {
        // remove leading delims and get first non-leading delim
        while ((index = to_split.find_first_of(delims)) == 0) {
            to_split = to_split.substr(1);
        }
        if (index == to_split.npos) {
            result.push_back(to_split);
            return result;
        } else {
            result.push_back(to_split.substr(0, index));
            to_split = to_split.substr(index + 1);
        }
    }
    return result;
}


void sort_controls(quantum_op_t *gate)
{
    std::vector<int> controls;
    for (int j = 0; j < 3; j++) {
        if (gate->ctrls[j] != -1) {
            controls.push_back(gate->ctrls[j]);
        }
    }
    sort(controls.begin(), controls.end());
    for (int j = 0; j < controls.size(); j++) {
        gate->ctrls[j] = controls[j];
    }
}


void sort_targets(quantum_op_t *gate)
{
    std::string name = std::string(gate->name);
    // only sort targets if gate is symmetric (i.e. U(t1,t2) = U(t2,t1))
    if (name == "swap" || name == "rzz" || name == "rxx" || name == "cswap") {
        int t1 = std::min(gate->targets[0], gate->targets[1]);
        int t2 = std::max(gate->targets[0], gate->targets[1]);
        gate->targets[0] = t1;
        gate->targets[1] = t2;
    }
}


void check_measurements(quantum_circuit_t* circuit)
{
    quantum_op_t* head = circuit->operations;
    bool seen_measurement = false;
    while(head != NULL) {
        if (head->type == op_measurement) {
            seen_measurement = true;
        }
        if (head->type == op_gate && seen_measurement) {
            circuit->has_intermediate_measurements = true;
            break;
        }
        head = head->next;
    }
    // TODO: also check if all qubits are measured
}


class QASMParser {

    private:

        // Gates specified in qelib1.inc
        // https://github.com/Qiskit/qiskit-terra/blob/main/qiskit/qasm/libs/qelib1.inc

        std::vector<std::string> qelib1_gates = 
        {
            "u3", "u2", "u1", "id", "u0", "u", "p", "x", "y", "z", "h", "s", "sdg",
            "t", "tdg", "rx", "ry", "rz", "sx", "sxdg",                                 // single qubit gates
            "cx", "cz", "cy", "swap", "ch", "crx", "cry", "crz", "cu1", "cp", "cu3", 
            "csx", "cu", "rxx", "rzz",                                                  // 2 qubit gates
            "ccx", "cswap", "rccx",                                                     // 3 qubit gates
            "rc3x", "c3x", "c3sx", "c3sqrtx",                                           // 4 qubit gates
            "c4x"                                                                       // 5 qubit gates
        }; 

        enum ins_type {
            comment, version_def, include, qreg, creg, barrier, measure, gate,
            gate_def, classical_cond, invalid
        };

        unsigned int current_line;

    public:

        typedef std::vector<std::pair<std::string, unsigned int>> registers_t;
    
        registers_t qregisters;
        registers_t cregisters;
    
        quantum_circuit_t *circuit;
        quantum_op_t *first_op;
        quantum_op_t *last_op;

        quantum_circuit_t* parse(char *filepath)
        {
            std::string path = std::string(filepath);
            //std::cout << "filepath = " << filepath << std::endl;
            std::string filename = path.substr(path.find_last_of("/\\") + 1);
            std::string circname = filename.substr(0, filename.find_last_of("."));
            std::ifstream infile(filepath);
            if (!infile.is_open()) {
                std::cerr << "Error opening file " << filepath << std::endl;
                exit(EXIT_FAILURE);
            }

            if (!infile.is_open()) {
                parse_error("Could not open the file!");
            }

            // create (blank) quantum circuit
            circuit = (quantum_circuit_t *) calloc(1, sizeof(quantum_circuit_t));
            first_op = (quantum_op_t *) calloc(1, sizeof(quantum_op_t));
            first_op->type = op_blank;
            first_op->next = NULL;
            circuit->operations = first_op;
            circuit->reversed_qubit_order = false;
            last_op = first_op;
            strcpy(circuit->name, circname.c_str());

            std::string line;
            current_line = 0;
            while (std::getline(infile, line))
            {
                parse_line(line);
                //std::cout << line << std::endl;
            }

            // if circuit has no intermediate measurements, make sure creg is
            // large enough to measure all qubits with measure_all
            check_measurements(circuit);
            if (!circuit->has_intermediate_measurements) {
                circuit->creg_size = std::max(circuit->creg_size, circuit->qreg_size);
            }
            circuit->creg = (bool *) calloc(circuit->creg_size, sizeof(bool));
            return circuit;
        }


        void parse_line(std::string line)
        {
            current_line++;
            if (line.size() == 0) {
                return;
            }
            auto args = split(line, " \t(");
            if (args.size() == 0) {
                return;
            }
            ins_type instruction_type;
            

            // 1. Get type of instruction
            instruction_type = get_instruction_type(&args[0]);

            // 2. Get arguments based on type of instruction
            switch (instruction_type)
            {
            case ins_type::comment:
                return;
            case ins_type::version_def:
                parse_version(line);
                return;
            case ins_type::include:
                // Includes are ignored for now
                return;
            case ins_type::barrier:
                // Currently, gates are executed in order they appear in the QASM file,
                // so barriers can be ignored.
                return;
            case ins_type::qreg:
                parse_reg(line, "qreg");
                break;
            case ins_type::creg:
                parse_reg(line, "creg");
                break;
            case ins_type::gate:
                parse_gate(args[0], line);
                break;
            case ins_type::measure:
                parse_measurement(line);
                break;
            case ins_type::gate_def:
                parse_error("Defining custom gates currently unsupported");
            case ins_type::classical_cond:
                parse_error("Classical conditioning currently unsupported");
            default:
                parse_error("Unsupported instruction type (" + line + ")");
            }
        }


        ins_type get_instruction_type(std::string *token)
        {
            if (token->rfind("//") == 0) {
                return ins_type::comment;;
            }
            else if (*token == "OPENQASM"){
                return ins_type::version_def;
            }
            else if (*token == "include") {
                return ins_type::include;
            }
            else if (*token == "qreg") {
                return ins_type::qreg;
            }
            else if (*token == "creg") {
                return ins_type::creg;
            }
            else if (*token == "measure") {
                return ins_type::measure;
            }
            else if (*token == "barrier") {
                return ins_type::barrier;
            }
            else if (*token == "gate") {
                return ins_type::gate_def;
            }
            else if (is_gate(token)) {
                return ins_type::gate;
            }
            else if (token->substr(0,2) == "if") {
                return ins_type::classical_cond;
            }
            else {
                return ins_type::invalid;
            }
        }


        void parse_version(std::string line)
        {
            auto args = split(line, " \t");

            if (args.size() > 0) {

                if (!args[0].compare("OPENQASM") || !args[1].substr(0,4).compare("2.0;")) {
                    return;
                }
                parse_error("Expected 'OPENQASM 2.0;', got something else instead"); // '" + args[1] + "'
            }
        }


        void parse_reg(std::string line, std::string type)
        {
            auto args = split(line, " []");
            if (args[0] != type) {
                parse_error("Expected '" + type + "', got '" + args[0] + "' instead");
            }
            if (args.size() < 3) {
                parse_error("Parsing error: Expected reg_name[size]");
            }

            if (type == "qreg") {
                qregisters.push_back({args[1], stoi(args[2])});
                circuit->qreg_size += stoi(args[2]);
            }
            else if (type == "creg") {
                cregisters.push_back({args[1], stoi(args[2])});
                circuit->creg_size += stoi(args[2]);
            }
            else {
                parse_error("Unexpected register type '" + type + "'");
            }
        }


        void parse_gate(std::string name, std::string line)
        {
            auto args = split(line, " ,[]()");
            if (args.size() == 0) {
                parse_error("Could not parse line");
            }

            // put measurement info into new quantum_op_t and append to circuit
            quantum_op_t* op = (quantum_op_t*) calloc(1, sizeof(quantum_op_t));
            op->type = op_gate;
            op->next = NULL;
            op->targets[0] = op->targets[1] = -1;
            op->ctrls[0] = op->ctrls [1] = op->ctrls[2] = -1;
            last_op->next = op;
            last_op = op;

            // single-qubit gates with no additional parameters
            if (name == "id" || name == "x" || name == "y" || name == "z" || 
                name == "h" || name == "s" || name == "sdg" || name == "t" || 
                name == "tdg" || name == "sx" || name == "sxdg") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // singe-qubit gates with a single angle
            else if (name == "rx" || name == "ry" || name == "rz" || name == "u0" ||
                     name == "u1" || name == "p") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[2], stoi(args[3]));
                    op->angle[0] = eval_math_expression(args[1]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // singe-qubit gates with two angles
            else if (name == "u2") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->angle[0] = eval_math_expression(args[1]);
                    op->angle[1] = eval_math_expression(args[2]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // single-qubit gates with three angles
            else if (name == "u3" || name == "u") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[4], stoi(args[5]));
                    op->angle[0] = eval_math_expression(args[1]);
                    op->angle[1] = eval_math_expression(args[2]);
                    op->angle[2] = eval_math_expression(args[3]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // two-qubit controlled gates with no additional parameters
            else if (name == "cx" || name == "cy" || name == "cz" || name == "ch" ||
                     name == "csx") {
                try {
                    strcpy(op->name, name.c_str());
                    op->targets[0] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // other two-qubit gates without additional parameters
            else if (name == "swap") {
                try {
                    strcpy(op->name, name.c_str());
                    op->targets[0] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->targets[1] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    sort_targets(op);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // two-qubit control gates with single angle
            else if (name == "crx" || name == "cry" || name == "crz" || name == "cp" ||
                     name == "cu1") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[4], stoi(args[5]));
                    op->ctrls[0] = get_seq_index(qregisters, args[2], stoi(args[3]));
                    op->angle[0] = eval_math_expression(args[1]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // two-qubit controlled gates with three angles
            else if (name == "cu3" || name == "cu") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[6], stoi(args[7]));
                    op->ctrls[0] = get_seq_index(qregisters, args[4], stoi(args[5]));
                    op->angle[0] = eval_math_expression(args[1]);
                    op->angle[1] = eval_math_expression(args[2]);
                    op->angle[2] = eval_math_expression(args[3]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // other two-qubit gate with two targets and a single angle
            else if (name == "rzz" || name == "rxx") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[4], stoi(args[5]));
                    op->targets[1] = get_seq_index(qregisters, args[2], stoi(args[3]));
                    op->angle[0] = eval_math_expression(args[1]);
                    sort_targets(op);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // three-qubit controlled gates with no additional parameters
            else if (name == "ccx" || name == "rccx") {
                try {
                    strcpy(op->name, name.c_str());
                    op->targets[0] = get_seq_index(qregisters, args[5], stoi(args[6]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    op->ctrls[1] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    sort_controls(op);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // other three-qubit gates with no additional parameters
            else if (name == "cswap") {
                try {
                    strcpy(op->name, name.c_str());
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    op->targets[0] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->targets[1] = get_seq_index(qregisters, args[5], stoi(args[6]));
                    sort_targets(op);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // four-qubit controlled gates with no additional parameters
            else if (name == "c3x" || name == "c3sx" || name == "c3sqrtx") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->targets[0] = get_seq_index(qregisters, args[7], stoi(args[8]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    op->ctrls[1] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->ctrls[2] = get_seq_index(qregisters, args[5], stoi(args[6]));
                    sort_controls(op);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            else {
                std::cout << current_line << ": " << line << std::endl;
                for (auto arg : args) {
                    std::cout << arg << ", ";
                }
                std::cout << std::endl;
                parse_error("Gate '" + name + "' currently unsupported");
            }
        }


        void parse_measurement(std::string line)
        {
            auto args = split(line, " []->");
            if (args.size() < 4) {
                parse_error("Expected more arguments for 'measure()'");
            }
            
            // put measurement info into new quantum_op_t
            quantum_op_t* op = (quantum_op_t*) malloc(sizeof(quantum_op_t));
            op->type = op_measurement;
            strcpy(op->name, "measure");
            op->targets[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
            op->meas_dest = -1;
            if (args.size() >= 4) {
                op->meas_dest = get_seq_index(cregisters, args[3], stoi(args[4]));
            }
            op->next = NULL;
            
            // append measurement to circuit
            last_op->next = op;
            last_op = op;
        }


        unsigned int get_seq_index(registers_t regs, std::string reg_name, unsigned int index)
        {
            unsigned int offset = 0;
            for (auto reg : regs) {
                if (reg.first == reg_name) {
                    return index + offset;
                }
                offset += reg.second;
            }
            parse_error("Register '" + reg_name + "' undefined");
            exit(EXIT_FAILURE);
        }


        bool is_gate(std::string *gate)
        {
            std::transform(gate->begin(), gate->end(), gate->begin(), tolower);
            
            return std::find(qelib1_gates.begin(), 
                             qelib1_gates.end(), 
                             *gate) 
                    != qelib1_gates.end();
        }


        std::string canonical_gate_name(std::string name)
        {
            if (name == "u0") return "id";
            if (name == "u1") return "p";
            if (name == "cu1") return "cp";
            if (name == "u3") return "u";
            if (name == "cu3") return "cu";
            if (name == "c3sqrtx") return "c3sx";
            return name;
        }


        void parse_error(std::string error) {
            std::cout << "Error in qasm file, line " << current_line << ": ";
            std::cout << error << std::endl;
            std::cout << "Parsing stopped" << std::endl;
            free_quantum_circuit(circuit);
            exit(EXIT_FAILURE);
        }

}; // QASMParser


quantum_circuit_t* parse_qasm_file(char *filepath)
{
    QASMParser parser;
    return parser.parse(filepath);
}


void print_quantum_op(quantum_op_t* op)
{
    if (op->type == op_measurement) {
        printf("measure(q%d)", op->targets[0]); 
        if (op->meas_dest >= 0) { 
            printf(" -> c%d", op->meas_dest);
        }
    }
    else if (op->type == op_gate) {
        printf("%s", op->name);
        for (int i = 0; i < 3; i++) {
            if (op->angle[i] != 0.0) {
                printf("_%lf", op->angle[i]);
            }
        }
        printf("(t={%d",op->targets[0]);
        if (op->targets[1] != -1) {
            printf(",%d",op->targets[1]);
        }
        printf("}");
        if (op->ctrls[0] != -1) {
            printf(",c={%d,%d,%d}", op->ctrls[0], op->ctrls[1], op->ctrls[2]);
        }
        printf(")");
    }
}


quantum_op_t* _get_swap_op(int q1, int q2)
{
    quantum_op_t* op = (quantum_op_t*) calloc(1, sizeof(quantum_op_t));
    op->type = op_gate;
    op->next = NULL;
    op->ctrls[0] = op->ctrls [1] = op->ctrls[2] = -1;
    strcpy(op->name, "swap");
    op->targets[0] = q1;
    op->targets[1] = q2;
    sort_targets(op);
    return op;
}


void insert_required_swaps(quantum_circuit_t *circuit)
{
    // For any controlled gate CU(c,t) with c > t (i.e. c below t in the DD), 
    // replace CU(c,t) with SWAP(c,t)CU(t,c)SWAP(c,t)

    quantum_op_t* head = circuit->operations->next;
    quantum_op_t* prev = circuit->operations; // first (blank) op always exists
    while(head != NULL) {
        if (head->type != op_gate) {
            head = head->next;
            continue;
        }

        // Check if target is below some controls
        bool should_swap = false;
        int c_index = 0;
        for (int j = 0; j < 3; j++) {
            if (head->ctrls[j] == -1) {
                break;
            }
            if (head->ctrls[j] > head->targets[0]) {
                should_swap = true;
            }
            if (head->ctrls[j] > head->ctrls[c_index]) {
                c_index = j;
            }
        }

        // If target is not below some controls, swap target with lowest control
        if (should_swap) {
            // swap control and target
            int tmp = head->targets[0];
            head->targets[0] = head->ctrls[c_index];
            head->ctrls[c_index] = tmp;

            // insert: prev -> swap1(t,c) -> head -> swap2(t,c) -> next
            quantum_op_t *next = head->next;
            quantum_op_t *swp1 = _get_swap_op(head->targets[0], head->ctrls[c_index]);
            quantum_op_t *swp2 = _get_swap_op(head->targets[0], head->ctrls[c_index]);
            sort_controls(head);
            sort_targets(head);
            prev->next = swp1;
            swp1->next = head;
            head->next = swp2;
            swp2->next = next;

            prev = swp2;
            head = next;
        } else {
            prev = head;
            head = head->next;
        }
    }
}


void order_cphase_gates(quantum_circuit_t *circuit)
{
    // Since for any controlled phase gate (cp, cz=cp(pi)) cp(c,t) = c(t,c),
    // change the order such that c < t
    // Note that crz(theta) is not symmetric

    quantum_op_t* head = circuit->operations;
    while (head != NULL) {
        if (head->type == op_gate) {
            std::string name = std::string(head->name);
            if (name == "cz" || name == "cp") {
                if (head->ctrls[0] > head->targets[0]) {
                    int tmp = head->targets[0];
                    head->targets[0] = head->ctrls[0];
                    head->ctrls[0] = tmp;
                }
            }
        }
        head = head->next;
    }
}


void reverse_order(quantum_circuit_t *circuit)
{
    // remap qubit index i to qreg_size-1-i
    
    quantum_op_t* head = circuit->operations;
    while (head != NULL) {
        if (head->type == op_gate || head->type == op_measurement) {
            for (int j = 0; j < 2; j++) {
                if (head->targets[j] != -1) {
                    head->targets[j] = (circuit->qreg_size - 1) - head->targets[j];
                }
            }
            for (int j = 0; j < 3; j++) {
                if (head->ctrls[j] != -1) {
                    head->ctrls[j] = (circuit->qreg_size - 1) - head->ctrls[j];
                }
            }
        }
        sort_controls(head);
        sort_targets(head);
        head = head->next;
    }
    circuit->reversed_qubit_order = true;
}


void optimize_qubit_order(quantum_circuit_t *circuit, bool allow_swaps)
{
    order_cphase_gates(circuit); // order phase gates before counting
    quantum_op_t* head = circuit->operations;
    int ctrls_below_target = 0;
    int ctrls_above_target = 0;
    while (head != NULL) {
        if (head->type == op_gate && head->ctrls[0] != -1) {
            if (head->ctrls[0] > head->targets[0]) {
                ctrls_below_target += 1;
            } else {
                ctrls_above_target += 1;
            }
        }
        head = head->next;
    }
    if (ctrls_below_target > ctrls_above_target) {
        reverse_order(circuit);
        order_cphase_gates(circuit); // order phase gates again
    }
    if (allow_swaps) {
        insert_required_swaps(circuit);
    }
}

quantum_op_t** circuit_as_array(quantum_circuit_t *circuit, bool return_non_empty, int *length)
{
    // loop over circuit to get lenght
    *length = 0;
    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        switch (op->type) {
            case op_gate:
            case op_measurement:
                *length += 1;
                break;
            default:
                break;
        }
        op = op->next;
    }

    if (*length == 0 && return_non_empty) {
        quantum_op_t **ops_array = (quantum_op_t**)malloc(sizeof(quantum_op_t*));
        quantum_op_t *id0 = (quantum_op_t*)malloc(sizeof(quantum_op_t));
        strcpy(id0->name, "id");
        id0->targets[0] = 0;
        ops_array[0] = id0;
        *length = 1;
        return ops_array;
    }
    
    // loop over circuit and put into array
    quantum_op_t **ops_array = (quantum_op_t**)malloc(sizeof(quantum_op_t*) * (*length));
    int i = 0;
    op = circuit->operations;
    while (op != NULL) {
        switch (op->type) {
            case op_gate:
            case op_measurement:
                ops_array[i] = op;
                i++;
                break;
            default:
                break;
        }
        op = op->next;
    }

    return ops_array;
}


void print_quantum_circuit(quantum_circuit_t *circuit)
{
    printf("qreg size: %d\n", circuit->qreg_size);
    printf("creg size: %d\n", circuit->creg_size);
    quantum_op_t *op = circuit->operations;
    while (op != NULL) {
        print_quantum_op(op);
        printf("\n");
        op = op->next;
    }
}

void fprint_creg(FILE *stream, quantum_circuit_t *circuit)
{
    // print in big-endian
    for (int i = circuit->creg_size-1; i >= 0 ; i--) {
        fprintf(stream, "%d", circuit->creg[i]);
    }
}


void free_quantum_circuit(quantum_circuit_t* circuit)
{
    quantum_op_t* head = circuit->operations;
    quantum_op_t* tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
    free(circuit->creg);
    free(circuit);
}
