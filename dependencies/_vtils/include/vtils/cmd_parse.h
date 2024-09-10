/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#ifndef VTILS_CMD_PARSE_h
#define VTILS_CMD_PARSE_h

#include "vtils/type_conv.h"

#include <iostream>
#include <map>
#include <set> 
#include <string>

#include <stdint.h>

namespace vtils {

class CmdParser {
public:
    CmdParser(int argc, char* argv[], int from=0);

    bool    option_set(std::string, bool exit_on_fail=false) const;
    bool    get(std::string, std::string&, bool exit_on_fail=false) const;

    template <class T> bool get(std::string, T&, bool exit_on_fail=false) const;
    
    void print_all_set_options(std::ostream&) const;

    std::string help;
private:
    std::set<std::string>               option_pool;
    std::map<std::string, std::string>  option_to_arg;
};

}   // vtils

#include "inl/cmd_parse.inl"

#endif  // VTILS_CMD_PARSE_h
