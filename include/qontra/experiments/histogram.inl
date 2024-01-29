/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "qontra/experiments.h"
#include "qontra/experiments/stats.h"

#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

#include <mpi.h>
#include <strings.h>

namespace qontra {

inline std::string
histogram_key_to_hex(histogram_key_t key) {
    std::ostringstream ss;
    for (size_t i = 0; i < key.size(); i++) {
        if (i) ss << " ";
        // For readability, every 2 bytes, we will insert a space.
        for (size_t j = 0; j < 4; j++) {
            uint64_t subword = key[i] >> (16*j);
            if (j) {
                // Check if the remaining bits are all 0.
                if (i == key.size()-1 && subword == 0) break;
                ss << " ";
            }
            uint16_t x = static_cast<uint16_t>(subword) & std::numeric_limits<uint16_t>::max();
            ss << std::hex << x;
        }
    }
    return ss.str();
}

template <class T> histogram_t<T>
histogram_reduce(const histogram_t<T>& hist) {
    if (!G_USE_MPI) return hist;
    int world_size, world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<std::pair<histogram_key_t, T>> hist_entries(hist.begin(), hist.end());

    histogram_t<T> reduced_hist;
    for (int r = 0; r < world_size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        // Broadcast the number of entries in hist.
        uint32_t n_entries = static_cast<uint32_t>(hist.size());
        MPI_Bcast(&n_entries, 1, MPI_UNSIGNED, r, MPI_COMM_WORLD);
        // Now send each entry:
        for (uint32_t i = 0; i < n_entries; i++) {
            histogram_key_t key;
            T value;
            if (world_rank == r) {
                auto tmp = hist_entries[i];
                key = tmp.first;
                value = tmp.second;
            }
            // Broadcast the size of tmp.
            uint32_t key_width = static_cast<uint32_t>(key.size());
            MPI_Bcast(&key_width, 1, MPI_UNSIGNED, r, MPI_COMM_WORLD);
            if (world_rank != r) {
                key = histogram_key_t(key_width);
            }
            // Now broadcast key and value.
            MPI_Bcast(&key[0], key_width, MPI_UNSIGNED_LONG, r, MPI_COMM_WORLD);
            MPI_Bcast(&value, 1, get_mpi_type<T>(), r, MPI_COMM_WORLD);
            reduced_hist[key] += value;
        }
    }
    return reduced_hist;
}

template <class T> inline T
histogram_sum_values(const histogram_t<T>& hist) {
    T s = static_cast<T>(0);
    for (const auto& p : hist) s += p.second;
    return s;
}

template <class T> histogram_t<double>
histogram_normalize(const histogram_t<T>& hist) {
    double s = static_cast<double>(histogram_sum_values(hist));
    histogram_t<double> norm_hist;
    for (const auto& p : hist) {
        norm_hist[p.first] = static_cast<double>(p.second) / s;
    }
    return norm_hist;
}

template <class T> inline size_t
histogram_get_max_key_width(const histogram_t<T>& hist) {
    size_t max_width = 0;
    for (const auto& p : hist) {
        size_t n_words = p.first.size();
        // The width is determined by the highest order bit in the last word.
        const size_t lower_order_width = (n_words-1) << 6;
        const size_t last_word_width = flsll(static_cast<long long>(p.first.back()));
        size_t width = lower_order_width + last_word_width;
        max_width = std::max(width, max_width);
    }
    return max_width;
}

template <class T> inline std::ostream&
operator<<(std::ostream& out, const histogram_t<T>& hist) {
    const size_t max_width = histogram_get_max_key_width(hist); // In bits
    const size_t hex_key_width_ignoring_spaces = max_width >> 3;
    // Spaces are introduced every four hex characters:
    const size_t hex_key_width_with_spaces = (hex_key_width_ignoring_spaces*5) / 4;
    // Add 3 characters of padding.
    const size_t key_width_with_padding = hex_key_width_with_spaces + 3;

    for (const auto& p : hist) {
        out << std::setw(key_width_with_padding) << histogram_key_to_hex(p.first)
            << " | "
            << p.second 
            << std::endl;
    }
    return out;
}

template <> inline std::ostream&
operator<< <double>(std::ostream& out, const histogram_t<double>& hist) {
    const size_t max_width = histogram_get_max_key_width(hist); // In bits
    const size_t hex_key_width_ignoring_spaces = max_width >> 3;
    // Spaces are introduced every four hex characters:
    const size_t hex_key_width_with_spaces = (hex_key_width_ignoring_spaces*5) / 4;
    // Add 3 characters of padding.
    const size_t key_width_with_padding = hex_key_width_with_spaces + 3;

    for (const auto& p : hist) {
        out << std::setw(key_width_with_padding) << histogram_key_to_hex(p.first)
            << " | "
            << std::setprecision(6) << p.second
            << std::endl;
    }
    return out;
}

}   // qontra
