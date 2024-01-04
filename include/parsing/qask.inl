/*
 *  author: Suhas Vittal
 *  date:   3 January 2024
 * */

#include <sstream>

namespace qontra {

inline bool
matches(std::string text, token_type t) {
    RE2::FullMatch(text, get_regex(t));
}

inline
QasmParser::read_tokens_onto_stack(std::string text) {
    std::istringstream iss(text);
    read_tokens_onto_stack(iss);
}

}   // qontra
