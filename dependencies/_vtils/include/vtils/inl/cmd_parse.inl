/*
 *  author: Suhas Vittal
 *  date:   29 December 2023
 * */

namespace vtils {

template <class T> bool
CmdParser::get(std::string opt, T& out, bool e) const {
    std::string tmp;
    if (!get(opt, tmp, e)) return false;
    out = conv<T>(tmp);
    return true;
}

}   // vtils
