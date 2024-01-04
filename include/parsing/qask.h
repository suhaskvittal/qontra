/*
 *  author: Suhas Vittal
 *  date:   3 January 2024
 *
 *  Why did I write my own lexer+parser? Bison and Flex drove me insane.
 * */

#ifndef QASK_h
#define QASK_h

#include <re2.h>

#include <iostream>
#include <tuple>
#include <vector>

namespace qontra {

enum class token_type {
    // Keywords:
    KW_repeat,
    // Not keywords:
    T_identifier,
    T_i_literal,
    T_f_literal,
    T_separator,
    T_brace_open,
    T_brace_close,
    T_modifier,
    T_whitespace,
    T_comment
};

RE2&    get_regex(token_type);
bool    matches(std::string, token_type);

// Each token has the following entries:
//  (1) token_type
//  (2) value (std::string)
typedef std::tuple<token_type, std::string> Token;

class QaskParser {
public:
    QaskParser();

    void read_tokens_onto_stack(std::string);
    void read_tokens_onto_stack(std::istream&);
private:
    std::vector<Token> token_stack;
};

}   // qontra

#include "qask.inl"

#endif  // QASK_h
