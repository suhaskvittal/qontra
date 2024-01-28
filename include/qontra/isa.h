/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#ifndef QONTRA_ISA_h
#define QONTRA_ISA_h

#include <map>
#include <string>
#include <vector>

namespace qontra {

enum class INSTRUCTION_TYPE {
    NOP                 = -1,
    QUANTUM_G1Q         = 0,
    QUANTUM_G2Q         = 1,
    EVENT_OR_OBS        = 2,
    REGISTER_LOGICAL    = 3,
    DATA_MOVEMENT       = 4,
    MEMORY_SHIFT        = 5,
    MEMORY_OFFSET       = 6,
    PERMISSIVE_BR_TYPE1 = 7,
    PERMISSIVE_BR_TYPE2 = 8,
    CALL                = 9
};

struct isa_data_t {
    INSTRUCTION_TYPE    inst_type;
    uint64_t            flags;

    // flags is a catch-all. Accessing specific data from the flags can
    // be done via the functions:
    bool    is_simd_like(void) { return flags & 0x1; }
    bool    error_precedes_op(void) { return flags & 0x2; }
    bool    apply_x_error_instead_of_depolarizing(void) { return flags & 0x4; }
};

const std::map<std::string, isa_data_t>&    isa(void);
bool                                        isa_is_initialized(void);

// This function checks if the ISA is initialized. If not, then it is initialized via
// build_isa_from_file (see below).
void    test_and_init_isa(void);

// Note: this function is case insensitive.
const isa_data_t&   isa_get(std::string);

std::vector<std::string>    isa_list_instructions(void);

// The ISA is specified via a file. By default, this file is in the macro
// QONTRA_ISA_FILE and is loaded on the first call to isa_get.
//
// The file must be a CSV file with the following structure per line:
//  
//  <name>,<instruction type (int)>,<flags (int)>
//
// i.e.
//  h,0,1
//  mshift,6,0
//  measure,0,7
//  cx,1,1
//
// This is case insensitive.
void    build_isa_from_file(std::string);

}   // qontra

#include "isa.inl"

#endif  // QONTRA_ISA_h
