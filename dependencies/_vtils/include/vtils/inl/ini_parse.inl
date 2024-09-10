/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

namespace vtils {

inline bool
IniParser::get(std::string key, std::string& value) const {
    return get("__ANON__", key, value);
}

inline bool
IniParser::get(std::string section, std::string key, std::string& value) const {
    if (!ini_map.count(section) || !ini_map.at(section).count(key)) return false;
    value = ini_map.at(section).at(key);
    return true;
}

template <class T> inline bool
IniParser::get(std::string key, T& value) const {
    return get("__ANON__", key, value);
}

template <class T> inline bool
IniParser::get(std::string section, std::string key, T& value) const {
    std::string tmp;
    if (get(section, key, tmp)) {
        value = conv<T>(tmp);
        return true;
    } else {
        return false;
    }
}

inline const std::map<std::string, std::map<std::string, std::string>>&
IniParser::get_ini_map() {
    return ini_map;
}

}   // vtils
