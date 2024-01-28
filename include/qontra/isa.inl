/*
 *  author: Suhas Vittal
 *  date:   28 January 2024
 * */

#include <algorithm>
#include <iostream>

namespace qontra {

inline void
test_and_init() {
    if (!isa_is_initialized()) {
#ifndef QONTRA_ISA_FILE
        std::cerr << "Preprocessor macro QONTRA_ISA_FILE is not set, so the ISA cannot be loaded." << std::endl;
        exit(1);
#else
        build_isa_from_file(QONTRA_ISA_FILE);
#endif
    }
}

inline const isa_data_t&
isa_get(std::string name) {
    test_and_init();

    std::transform(name.begin(), name.end(), name.begin(),
            [] (char x) {
                return std::tolower(x);
            });
    if (!isa().count(name)) {
        std::cerr << "Instruction \"" << name << "\" is not in ISA!" << std::endl;
        exit(1);
    }
    return isa().at(name);
}

inline std::vector<std::string>
isa_list_instructions() {
    test_and_init();

    std::vector<std::string> arr;
    for (const auto& p : isa()) arr.push_back(p.first);
    return arr;
}


}   // qontra
