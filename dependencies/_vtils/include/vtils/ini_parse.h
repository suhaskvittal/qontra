/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#ifndef VTILS_INI_PARSE_h
#define VTILS_INI_PARSE_h

#include "vtils/type_conv.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

namespace vtils {

class IniParser {
public:
    IniParser(std::string file);

    bool get(std::string key, std::string& value) const;
    bool get(std::string section, std::string key, std::string& value) const;

    // This reads the string at [section][key] and converts to a type T using
    // the FUNC.
    template <class T> bool             get(std::string key, T& value) const;
    template <class T> bool             get(std::string section, std::string key, T& value) const;

    const std::map<std::string, std::map<std::string, std::string>>&
        get_ini_map(void);
private:
    void parse(std::ifstream&);

    std::map<std::string, std::map<std::string, std::string>> ini_map;
};

}   // vtils

#include "inl/ini_parse.inl"

#endif  // VTILS_INI_PARSE_h
