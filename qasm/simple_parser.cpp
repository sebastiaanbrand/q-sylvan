#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
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
        index = to_split.find_first_of(delims);
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
        std::vector<std::pair<std::string, unsigned int>> qregisters;
        std::vector<std::pair<std::string, unsigned int>> cregisters;


        gate_info_t* parse(char *filepath)
        {
            std::cout << "Parsing file " << std::string(filepath) << std::endl;

            std::ifstream infile(filepath);

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

            return NULL;
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
                std::cout << current_line << ": " << line << std::endl;
                parse_measurement(line);
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
            // TODO
        }


        unsigned int get_sequential_index(std::string reg, unsigned int index)
        {
            // TODO
            return 0;
        }

}; // QASMParser


gate_info_t* parse_qasm_file(char *filepath)
{
    QASMParser parser;
    return parser.parse(filepath);
}
