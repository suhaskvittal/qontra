/*
 *  author: Suhas Vittal
 *  date:   30 May 2024
 * */

#ifndef PLACC_CX_h
#define PLACC_CX_h

#include <qes.h>

#include <set>
#include <vector>

namespace placc {

struct CXLayer {
    bool test_and_add(uint64_t, uint64_t);
    bool test_and_add(int64_t, int64_t);
    // int64_t args to conform to qes.
    std::set<int64_t> in_use;
    std::vector<int64_t> operands;
};

class CXManager {
public:
    CXManager(void);

    size_t get_depth(void);
    void push_back_cx(uint64_t, uint64_t);
    void flush(qes::Program<>&);
    qes::Program<> flush(void);
private:
    std::vector<CXLayer> layers;
    std::map<uint64_t, size_t> operand_layer_map;
};

}   // placc

#include "inl/cx.inl"

#endif  // PLACC_CX_h
