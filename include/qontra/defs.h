/* 
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DEFS_h
#define DEFS_h

#include <memory>
#include <vector>
#include <stdint.h>

typedef double fp_t;

typedef std::vector<fp_t>   poly_t;

template <class T> using sptr=std::shared_ptr<T>;
template <class T> using uptr=std::unique_ptr<T>;

#endif
