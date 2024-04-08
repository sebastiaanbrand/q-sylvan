#include <string>

#include "eval_expr.hpp"
#include "exprtk.hpp"

typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> math_expr_t;
typedef exprtk::parser<double> math_parser_t;

double eval_math_expression(std::string expression_string)
{
    symbol_table_t symbol_table;
    symbol_table.add_constants(); // this add pi and other constants

    math_expr_t expression;
    expression.register_symbol_table(symbol_table);

    math_parser_t math_parser;
    math_parser.compile(expression_string, expression);

    return expression.value();
}