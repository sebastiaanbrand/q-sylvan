#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#include "simple_parser.h"
#include "parse_math/eval_expr.hpp"

/**
 * Minor shortcomings of this parser:
 * 1) Assumes only one instruction per line
 * 2) Ignores includes
 * 3) Assumes that for given rotation angles (e.g. rz(pi/4)), the expression 
 *    between the brackets doesn't contain additional brackets.
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


class QASMParser {

    private:
        // Gates specified in qelib1.inc
        // https://github.com/Qiskit/qiskit-terra/blob/main/qiskit/qasm/libs/qelib1.inc
        std::vector<std::string> qelib1_gates = 
        {"u3", "u2", "u1", "id", "u0", "u", "p", "x", "y", "z", "h", "s", "sdg",
        "t", "tdg", "rx", "ry", "rz", "sx", "sxdg", // single qubit gates
        "cx", "cz", "cy", "swap", "ch", "crx", "cry", "crz", "cu1", "cp", "cu3", 
        "csx", "cu", "rxx", "rzz", // 2 qubit gates
        "ccx", "cswap", "rccx", // 3 qubit gates
        "rc3x", "c3x", "c3sx", "c3sqrtx", // 4 qubit gates
        "c4x"}; // 5 qubit gates

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
            std::cout << "Parsing file " << std::string(filepath) << std::endl;
            std::ifstream infile(filepath);

            // create (blank) quantum circuit
            circuit = (quantum_circuit_t *) malloc(sizeof(quantum_circuit_t));
            first_op = (quantum_op_t *) malloc(sizeof(quantum_op_t));
            first_op->type = op_blank;
            first_op->next = NULL;
            circuit->operations = first_op;
            last_op = first_op;

            std::string line;
            current_line = 0;
            while (std::getline(infile, line))
            {
                parse_line(line);
            }

            // remove (blank) first operation
            circuit->operations = circuit->operations->next;
            free(first_op);

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
                parse_error("Unsupported instruction type");
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
                if (args[1] == "2.0;") {
                    return;
                }
                std::cerr << "WARNING: expected OPENQASM version 2.0, "
                      << "got \"" << args[1] << "\" instead" << std::endl;
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
            }
            else if (type == "creg") {
                cregisters.push_back({args[1], stoi(args[2])});
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
            op->ctrls[0] = op->ctrls [1] = op->ctrls[2] = -1;
            last_op->next = op;
            last_op = op;

            // single-qubit gates with no additional parameters
            if (name == "id" || name == "x" || name == "y" || name == "z" || 
                name == "h" || name == "s" || name == "sdg" || name == "t" || 
                name == "tdg" || name == "sx" || name == "sxdg") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->target = get_seq_index(qregisters, args[1], stoi(args[2]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // singe-qubit gates with a single angle
            else if (name == "rx" || name == "ry" || name == "rz" || name == "u0" ||
                     name == "u1" || name == "p") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->target = get_seq_index(qregisters, args[2], stoi(args[3]));
                    op->angle[0] = eval_math_expression(args[1]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // singe-qubit gates with two angles
            else if (name == "u2") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->target = get_seq_index(qregisters, args[3], stoi(args[4]));
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
                    op->target = get_seq_index(qregisters, args[4], stoi(args[5]));
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
                    op->target = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // other two-qubit gates without additional parameters
            else if (name == "swap") {
                try {
                    strcpy(op->name, name.c_str());
                    op->target = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // two-qubit controlled gates with single angle
            else if (name == "crx" || name == "cry" || name == "crz" || name == "cp" ||
                     name == "cu1") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->target = get_seq_index(qregisters, args[4], stoi(args[5]));
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
                    op->target = get_seq_index(qregisters, args[6], stoi(args[7]));
                    op->ctrls[0] = get_seq_index(qregisters, args[4], stoi(args[5]));
                    op->angle[0] = eval_math_expression(args[1]);
                    op->angle[1] = eval_math_expression(args[2]);
                    op->angle[2] = eval_math_expression(args[3]);
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // three-qubit controlled gates with no additional parameters
            else if (name == "ccx") {
                try {
                    strcpy(op->name, name.c_str());
                    op->target = get_seq_index(qregisters, args[5], stoi(args[6]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    op->ctrls[1] = get_seq_index(qregisters, args[3], stoi(args[4]));
                } catch (...) {
                    parse_error("Error parsing arguments of gate " + name);
                }
            }
            // four-qubit controlled gates with no additional parameters
            else if (name == "c3x" || name == "c3sx" || name == "c3sqrtx") {
                try {
                    strcpy(op->name, canonical_gate_name(name).c_str());
                    op->target = get_seq_index(qregisters, args[7], stoi(args[8]));
                    op->ctrls[0] = get_seq_index(qregisters, args[1], stoi(args[2]));
                    op->ctrls[1] = get_seq_index(qregisters, args[3], stoi(args[4]));
                    op->ctrls[2] = get_seq_index(qregisters, args[5], stoi(args[6]));
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
            auto args = split(line, " []");
            if (args.size() < 4) {
                parse_error("Expected more arguments to 'measure'");
            }
            
            // put measurement info into new quantum_op_t
            quantum_op_t* op = (quantum_op_t*) malloc(sizeof(quantum_op_t));
            op->type = op_measurement;
            strcpy(op->name, "measure");
            op->target = get_seq_index(qregisters, args[1], stoi(args[2]));
            op->meas_dest = -1;
            if (args.size() >= 5) {
                op->meas_dest = get_seq_index(cregisters, args[4], stoi(args[5]));
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
            std::cout << "Parsing error on line " << current_line << ": ";
            std::cout << error << std::endl;
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
        printf("measure(q%d)", op->target); 
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
        printf("(%d", op->target);
        if (op->ctrls[0] >= 0) {
            printf(",c=%d,%d,%d)", op->ctrls[0], op->ctrls[1], op->ctrls[2]);
        } else {
            printf(")");
        }
    }
}


void print_quantum_circuit(quantum_circuit_t* circuit)
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


void free_quantum_circuit(quantum_circuit_t* circuit)
{
    quantum_op_t* head = circuit->operations;
    quantum_op_t* tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
    free(circuit);
}
