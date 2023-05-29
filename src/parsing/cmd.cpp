/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#include "parsing/cmd.h"

namespace qontra {

const std::regex no_arg_opt("-([A-Za-z]+)");
const std::regex w_arg_opt("--([A-Za-z-]+)");

CmdParser::CmdParser(int argc, char* argv[])
    :option_pool(),
    option_to_arg()
{
    int ptr = 1;
    while (ptr < argc) {
        std::string s(argv[ptr]);
        std::smatch m;
        if (std::regex_match(s, m, no_arg_opt)) {
            option_pool.insert(m[1]);
        } else if (std::regex_match(s, m, w_arg_opt)) {
            option_pool.insert(m[1]);
            option_to_arg[m[1]] = argv[++ptr];
        }
        ptr++;
    }
}

}   // qontra
