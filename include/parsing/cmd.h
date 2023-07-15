/*
 *  author: Suhas Vittal
 *  date:   29 May 2023
 * */

#ifndef PARSING_CMD_h
#define PARSING_CMD_h

#include "defs.h"

#include <map>
#include <regex>
#include <set>
#include <string>
#include <vector>

#include <stdint.h>

namespace qontra {

class CmdParser {
public:
    CmdParser(int argc, char* argv[]);

    bool option_set(std::string opt) { return option_pool.count(opt); }

    bool get_string(std::string opt, std::string& out) {
        if (!option_set(opt))   return false;
        out = option_to_arg[opt]; 
        return true;
    }

    bool get_float(std::string opt, fp_t& out) {
        if (!option_set(opt))   return false;
        out = std::stof(option_to_arg[opt]); 
        return true;
    }

    bool get_int32(std::string opt, int32_t& out) {
        if (!option_set(opt))   return false;
        out = std::stoi(option_to_arg[opt]);
        return true;
    }

    bool get_uint32(std::string opt, uint32_t& out) {
        if (!option_set(opt))   return false;
        out = std::stoi(option_to_arg[opt]); 
        return true;
    }
    bool get_uint64(std::string opt, uint64_t& out) {
        if (!option_set(opt))   return false;
        out = std::stoll(option_to_arg[opt]); 
        return true;
    }
private:
    std::set<std::string>               option_pool;
    std::map<std::string, std::string>  option_to_arg;
};

}   // qontra

#endif  // PARSING_CMD_h
