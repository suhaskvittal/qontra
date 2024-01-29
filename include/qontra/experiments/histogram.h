/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#ifndef QONTRA_EXPERIMENTS_HISTOGRAM_h
#define QONTRA_EXPERIMENTS_HISTOGRAM_h

#include <map>
#include <vector>
#include <string>

#include <stdint.h>

namespace qontra {
    
typedef std::vector<uint64_t> histogram_key_t;

template <class NUMBER> using histogram_t=std::map<histogram_key_t, NUMBER>;

std::string histogram_key_to_hex(histogram_key_t);

template <class T> histogram_t<T>       histogram_reduce(const histogram_t<T>&);
template <class T> T                    histogram_sum_values(const histogram_t<T>&);
template <class T> histogram_t<double>  histogram_normalize(const histogram_t<T>&);
template <class T> size_t               histogram_get_max_key_width(const histogram_t<T>&);

template <class T>
std::ostream& operator<<(std::ostream&, const histogram_t<T>&);

}   // qontra

#include "histogram.inl"

#endif  // QONTRA_EXPERIMENTS_HISTOGRAM_h
