/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include "qontra/isa.h"

#include <fstream>
#include <iostream>

namespace qontra {

static bool ISA_INITIALIZED = false;
static std::map<std::string, isa_data_t> ISA;

const std::map<std::string, isa_data_t>& isa() {
    return ISA; 
}

bool isa_is_initialized() {
    return ISA_INITIALIZED;
}

void
build_isa_from_file(std::string filename) {
    std::ifstream fin(filename);

    std::string ln;
    while (std::getline(fin, ln)) {
        if (ln.empty()) continue;
        // Parse the line.
        size_t pptr, cptr;
        // Get instruction name:
        pptr = 0;
        cptr = ln.find(",");
        if (cptr == std::string::npos) {
            std::cerr << "Failed to parse instruction name in ISA." << std::endl;
            exit(1);
        }
        std::string inst_name = ln.substr(pptr, cptr-pptr);
        // Get instruction type.
        pptr = cptr+1;
        cptr = ln.find(",", pptr);
        if (cptr == std::string::npos) {
            std::cerr << "Failed to find instruction type when reading instruction \""
                << inst_name << "\"." << std::endl;
            exit(1);
        }
        INSTRUCTION_TYPE inst_type = static_cast<INSTRUCTION_TYPE>(std::stoi(ln.substr(pptr, cptr-pptr)));
        // Get instruction flags.
        pptr = cptr+1;
        uint64_t flags = std::stoull(ln.substr(pptr));

        isa_data_t dat = { inst_type, flags };
        ISA[inst_name] = dat;
    }
    ISA_INITIALIZED = true;
}

}   // qontra
