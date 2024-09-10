/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

namespace vtils {

template <class T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, bool>>
inline T conv(std::string s) {
    return static_cast<T>(stoull(s));
}

template <class T, std::enable_if_t<std::is_integral<T>::value && !std::is_unsigned<T>::value, bool>>
inline T conv(std::string s) {
    return static_cast<T>(stoll(s));
}

template <class T, std::enable_if_t<std::is_floating_point<T>::value, bool>>
inline T conv(std::string s) {
    return static_cast<T>(stold(s));
}

}   // vtils
