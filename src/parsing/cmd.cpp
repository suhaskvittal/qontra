/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#include "parsing/cmd.h"

#include <regex>

const std::regex no_arg_opt("-([1-9A-Za-z_-]+)");
const std::regex w_arg_opt("--([1-9A-Za-z_-]+)");

CmdParser::CmdParser(int argc, char* argv[])
    :option_pool(),
    option_to_arg()
{
    int ptr = 1;
    while (ptr < argc) {
        std::string s(argv[ptr]);
        std::smatch m;
        if (std::regex_match(s, m, w_arg_opt)) {
            option_pool.insert(m[1]);
            option_to_arg[m[1]] = argv[++ptr];
        } else if (std::regex_match(s, m, no_arg_opt)) {
            option_pool.insert(m[1]);
        }
        ptr++;
    }
}
