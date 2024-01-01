/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

inline bool
CmdParser::option_set(std::string opt) {
    return option_pool.count(opt); 
}

inline bool
CmdParser::get_string(std::string opt, std::string& out) {
    if (!option_set(opt))   return false;
    out = option_to_arg[opt]; 
    return true;
}

inline bool
CmdParser::get_float(std::string opt, double& out) {
    if (!option_set(opt))   return false;
    out = std::stof(option_to_arg[opt]); 
    return true;
}

inline bool
CmdParser::get_int32(std::string opt, int32_t& out) {
    if (!option_set(opt))   return false;
    out = std::stoi(option_to_arg[opt]);
    return true;
}

inline bool
CmdParser::get_int64(std::string opt, int64_t& out) {
    if (!option_set(opt))   return false;
    out = std::stoll(option_to_arg[opt]); 
    return true;
}

inline bool
CmdParser::get_uint32(std::string opt, uint32_t& out) {
    if (!option_set(opt))   return false;
    out = std::stoul(option_to_arg[opt]); 
    return true;
}
inline bool
CmdParser::get_uint64(std::string opt, uint64_t& out) {
    if (!option_set(opt))   return false;
    out = std::stoull(option_to_arg[opt]); 
    return true;
}

inline void
CmdParser::print_all_set_options(std::ostream& out) {
    for (std::string opt : option_pool) {
        out << opt;
        if (option_to_arg.count(opt)) {
            out << ": " << option_to_arg[opt];
        }
        out << "\n";
    }
}
