/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "vtils/ini_parse.h"

#include <regex>

namespace vtils {

IniParser::IniParser(std::string file)
    :ini_map()
{
    std::ifstream fin(file);
    parse(fin);
}

void
IniParser::parse(std::ifstream& fin) {
    // Read line by line
    // 
    // Rules:
    //  (1) ignore whitespace while reading
    //  (2) ; and # are comments
    //  (3) characters inside quotes considered raw characters.
    //  (4) delimited characters circumvent the above rules.
    std::string ln;

    std::string token;
    std::string curr_section = "__ANON__";
    std::string curr_key = "";
    while (std::getline(fin, ln)) {
        // Check if this matches a section header.
        bool in_comment = false;
        bool in_section_header = false;
        bool in_quote = false;
        bool in_value = false;
        bool escape_set = false;

        for (size_t i = 0; i < ln.size(); i++) {
            char c = ln[i];
            if (in_comment) continue;

            bool escape_set_on_this_char = false;
            if (!escape_set && !in_quote) {
                if (!in_section_header && (c == ' ' || c == '\t')) {
                    continue;
                } else if (c == ';' || c == '#') {
                    in_comment = true;
                } else if (c == '\\') {
                    escape_set_on_this_char = true;
                } else if (!in_section_header && c == '[') {
                    in_section_header = true;
                } else if (in_section_header && c == ']') {
                    in_section_header = false;
                    // Token becomes section header.
                    curr_section = std::move(token);
                    token.clear();
                } else if (c == '=') {
                    in_value = true;
                    // Current token is a key.
                    curr_key = std::move(token);
                    token.clear();
                } else if (c == '\"' || c == '\'') {
                    in_quote = true;
                } else {
                    token.push_back(c);
                }
            } else if (!escape_set && in_quote) {
                if (c == '\"' || c == '\'') {
                    in_quote = false;
                } else {
                    token.push_back(c);
                }
            } else {
                token.push_back(c);
            }
            escape_set = escape_set_on_this_char;
        }
        // Now that we are done with the line, handle anything remaining.
        if (!escape_set && !in_quote) {
            if (in_value) {
                ini_map[curr_section][curr_key] = std::move(token);
            } else if (in_comment) {
                in_comment = false;
            }
            token.clear();
        } else {
            token.push_back('\n');
        }
    }
}

}   // vtils
