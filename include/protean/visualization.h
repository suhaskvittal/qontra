/*
 *  author: Suhas Vittal
 *  date:   28 December 2023
 * */

#ifndef PROTEAN_VISUALIZATION_h
#define PROTEAN_VISUALIZATION_h

#include "protean/network.h"

#include <graphviz/gvc.h>

namespace qontra {
namespace protean {

void write_to_dot(PhysicalNetwork&);

}   // protean
}   // qontra

#endif  // PROTEAN_VISUALIZATION_h
