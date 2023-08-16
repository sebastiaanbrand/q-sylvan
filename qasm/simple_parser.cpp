#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>

#include "simple_parser.h"

/**
 * Minor shortcomings of this parser:
 * 1) Assumes only one instruction per line
 * 2) Ignores includes
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
        enum gate {
            I,X,Y,Z,H
        };

        enum ins_type {
            comment, version_def, include, qreg, creg, barrier, measure, gate
        };

        unsigned int current_line;

    public:
        typedef std::vector<std::pair<std::string, unsigned int>> registers_t;
        registers_t qregisters;
        registers_t cregisters;
        quantum_op_t* first_op;
        quantum_op_t* last_op;


        quantum_op_t* parse(char *filepath)
        {
            std::cout << "Parsing file " << std::string(filepath) << std::endl;
            std::ifstream infile(filepath);

            // create (blank) initial element for LL which represents the circuit
            first_op = (quantum_op_t *) malloc(sizeof(quantum_op_t));
            first_op->type = op_blank;
            first_op->next = NULL;
            last_op = first_op;

            std::string line;
            current_line = 0;
            while (std::getline(infile, line))
            {
                parse_line(line);
            }

            // temp summary (TODO: remove)
            std::cout << "Quantum registers:" << std::endl;
            for (auto qreg : qregisters) {
                std::cout << qreg.first << ", " << qreg.second << std::endl;
            }
            std::cout << "Classical registers:" << std::endl;
            for (auto creg : cregisters) {
                std::cout << creg.first << ", " << creg.second << std::endl;
            }
            std::cout << "Circuit:" << std::endl;
            print_quantum_ops(first_op);

            return first_op;
        }


        void parse_line(std::string line)
        {
            current_line++;
            std::stringstream line_ss(line);
            std::string token;
            ins_type instruction_type;
            

            // 1. Get type of instruction
            if (std::getline(line_ss, token, ' ')) {
                instruction_type = get_instruction_type(token);
            }

            // 2. Get arguments based on type of instruction
            switch (instruction_type)
            {
            case ins_type::comment:
                return;
            case ins_type::version_def:
                parse_version(line_ss);
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
                std::cout << current_line << ": " << line << std::endl;
                parse_gate(line);
                break;
            case ins_type::measure:
                parse_measurement(line);
                break;
            default:
                std::cerr << "Invalid instruction type: " << std::endl;
                std::cerr << current_line << ": " << line << std::endl;
                exit(EXIT_FAILURE);
            }
        }


        ins_type get_instruction_type(std::string token)
        {
            if (token.rfind("//") == 0) {
                return ins_type::comment;;
            }
            else if (token == "OPENQASM"){
                return ins_type::version_def;
            }
            else if (token == "include") {
                return ins_type::include;
            }
            else if (token == "qreg") {
                return ins_type::qreg;
            }
            else if (token == "creg") {
                return ins_type::creg;
            }
            else if (token == "measure") {
                return ins_type::measure;
            }
            else if (token == "barrier") {
                return ins_type::barrier;
            }
            else {
                return ins_type::gate;
            }
        }


        void parse_version(std::stringstream &args)
        {
            std::string token = "";
            if (std::getline(args, token, ' ')) {
                if (token == "2.0;") {
                    return;
                }
            }
            std::cerr << "WARNING: expected OPENQASM version 2.0, "
                    << "got \"" << token << "\" instead" << std::endl;
        }


        void parse_reg(std::string line, std::string type)
        {
            auto args = split(line, " []");
            if (args[0] != type) {
                std::cerr << "Parsing error: expected \"" << type << "\", got \""
                          << args[0] << "\" instead" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (args.size() < 3) {
                std::cerr << "Parsing error: Expected name[num] after " << type
                          << ", got \"" << line << "\" instead" << std::endl; 
            }

            if (type == "qreg") {
                qregisters.push_back({args[1], stoi(args[2])});
            }
            else if (type == "creg") {
                cregisters.push_back({args[1], stoi(args[2])});
            }
            else {
                std::cerr << "Error: Unexpected register type \"" << type << "\""
                          << std::endl;
                exit(EXIT_FAILURE);
            }
        }


        void parse_gate(std::string line)
        {
            // TODO
        }


        void parse_measurement(std::string line)
        {
            auto args = split(line, " []");
            if (args.size() < 7) {
                std::cerr << "Parsing error on line " << current_line << ": "
                          << "Expected more arguments to 'measure'" << std::endl;
                exit(EXIT_FAILURE);
            }
            
            // put measurement info into new quantum_op_t
            quantum_op_t* op = (quantum_op_t*) malloc(sizeof(quantum_op_t));
            op->type = op_measurement;
            strcpy(op->name, "measure");
            op->target = get_seq_index(qregisters, args[1], stoi(args[2]));
            op->meas_dest = get_seq_index(cregisters, args[4], stoi(args[5]));
            op->next = NULL;
            
            // append to circuit
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
            std::cerr << "Parsing error on line " << current_line << ": "
                      << "register '" << reg_name << "' undefined" << std::endl;
            exit(EXIT_FAILURE);
        }

}; // QASMParser


quantum_op_t* parse_qasm_file(char *filepath)
{
    QASMParser parser;
    return parser.parse(filepath);
}


void print_quantum_op(quantum_op_t* op)
{
    if (op->type == op_measurement) {
        printf("measure(q%d) -> c%d", op->target, op->meas_dest);
    }
    else if (op->type == op_gate) {
        printf("%s_%lf(t=%d, c=%d,%d,%d)", op->name, op->angle, op->target, 
                                           op->controls[0], op->controls[1], 
                                           op->controls[2]);
    }
}


void print_quantum_ops(quantum_op_t* head)
{
    while (head != NULL) {
        print_quantum_op(head);
        printf("\n");
        head = head->next;
    }
}


void free_quantum_ops(quantum_op_t* head)
{
    quantum_op_t* tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp);
    }
}
