/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#ifndef VTILS_TYPE_CONV_h
#define VTILS_TYPE_CONV_h

#include <string>
#include <type_traits>

namespace vtils {

template <class T, std::enable_if_t<std::is_integral<T>::value && std::is_unsigned<T>::value, bool> = true>
T conv(std::string);

template <class T, std::enable_if_t<std::is_integral<T>::value && !std::is_unsigned<T>::value, bool> = true>
T conv(std::string);

template <class T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
T conv(std::string);

}   // vtils

#include "inl/type_conv.inl"

#endif  // VTILS_TYPE_CONV_h
