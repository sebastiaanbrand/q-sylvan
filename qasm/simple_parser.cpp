#include <iostream>
#include <fstream>
#include <sstream>
#include "simple_parser.h"

/**
 * Minor shortcomings of this parser:
 * 1) Assumes only one instruction per line
 * 2) Ignores includes
*/

enum ins_type {
    comment,
    version_def,
    include,
    qreg,
    creg,
    barrier,
    measure,
    gate
};

enum gate {
    I,
    X,
    Y,
    Z,
    H
};



ins_type get_instruction_type(std::string token, int line_no)
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

void parse_qreg(std::stringstream &args)
{
    // TODO
}

void parse_creg(std::stringstream &args)
{
    // TODO
}

void parse_gate(std::string gate, std::stringstream &args)
{
    // TODO
}

void parse_measurement(std::stringstream &args)
{
    // TODO
}


void parse_line(std::string line, int line_no)
{
    std::cout << line_no << ": " << line << std::endl;

    std::stringstream line_ss(line);
    std::string token;
    ins_type instruction_type;

    // 1. Get type of instruction
    if (std::getline(line_ss, token, ' ')) {
        instruction_type = get_instruction_type(token, line_no);
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
        parse_qreg(line_ss);
        break;
    case ins_type::creg:
        parse_creg(line_ss);
        break;
    case ins_type::gate:
        parse_gate(token, line_ss);
        break;
    case ins_type::measure:
        parse_measurement(line_ss);
    default:
        std::cerr << "Invalid instruction type" << std::endl;
        exit(EXIT_FAILURE);
    }

}

gate_info_t* parse_qasm_file(char *filepath)
{
    std::cout << "Parsing file " << std::string(filepath) << std::endl;

    std::ifstream infile(filepath);

    std::string line;
    int line_no = 0;
    while (std::getline(infile, line))
    {
        parse_line(line, ++line_no);
    }

    return NULL;
}
