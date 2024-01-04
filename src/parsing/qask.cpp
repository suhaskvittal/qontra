/*
 *  author: Suhas Vittal
 *  date:   3 January 2024
 * */

#include "parsing/qask.h"

#include <map>
#include <set>

namespace qontra {

#define PAIR(x, y) std::make_pair(token_type::##x, RE2(y))

const std::map<token_type, RE2> REGEX_MAP{
    PAIR(KW_repeat, "repeat"),
    PAIR(T_identifier, "[A-Za-z]+"),
    PAIR(T_i_literal, R"_(\d+)_"),
    PAIR(T_f_literal, R"_(\d*.\d+|\d+.\d*)_"),
    PAIR(T_separator, ";"),
    PAIR(T_brace_open, "("),
    PAIR(T_brace_close, ")"),
    PAIR(T_modifier, "@"),
    PAIR(T_whitespace, R"_(\s+)_"),
    PAIR(T_comment, R"_(#.*?\n)_")
};

const std::vector<token_type> token_order{
    token_type::KW_repeat,
    token_type::T_identifier,
    token_type::T_i_literal,
    token_type::T_f_literal,
    token_type::T_separator,
    token_type::T_brace_open,
    token_type::T_brace_close,
    token_type::T_modifier,
    token_type::T_whitespace,
    token_type::T_comment
};

inline RE2&
get_regex(token_type t) {
    return REGEX_MAP[t];
}

QaskParser::QaskParser()
    :token_stack()
{}

void
QaskParser::read_tokens_onto_stack(std::istream& input) {
    std::string prev_symbol;
    std::string curr_symbol;
    token_type symbol_token = token_type::T_identifier;

    char c;
    while (!input.eof()) {
        prev_symbol = curr_symbol;
        // Update curr symbol.
        c = input.get();
        curr_symbol.push_back(c);
        // Recheck token regex.
        for (token_t t : token_order) {
            if (matches(curr_symbol, t)) {
                symbol_token = t;
                goto match_found;
            }
        }
        // If nothing matches the token, then push back prev_symbol
        // onto the token_stack.
        token_stack.push_back(std::make_tuple(symbol_token, prev_symbol));
        // Reset curr_symbol to just c.
        curr_symbol.clear();
        curr_symbol.push_back(c);
match_found:
        // Nothing to be done.
    }
    // Place the existing symbol onto the stack.
    token_stack.push_back(std::make_tuple(symbol_token, curr_symbol));
}

}   // qontra

