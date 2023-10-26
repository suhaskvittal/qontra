/*
 *  author: Suhas Vittal
 *  date:   26 October 2023
 * */

#ifndef PROTEAN_UTILS_h
#define PROTEAN_UTILS_h

#include "instruction.h"
#include "protean/representation.h"

#include <set>
#include <vector>

namespace qontra {
namespace protean {

schedule_t
write_memory_experiment(css_code_data_t, uint rounds, bool is_memory_x);

schedule_t
write_syndrome_extraction_ops(css_code_data_t);

}   // protean
}   // qontra

#endif  // PROTEAN_UTILS_h
