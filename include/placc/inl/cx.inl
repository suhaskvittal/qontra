/*
 *  author: Suhas Vittal
 *  date:   30 May 2024
 * */

namespace placc {

inline bool
CXLayer::test_and_add(uint64_t x, uint64_t y) {
    return test_and_add(static_cast<int64_t>(x), static_cast<int64_t>(y));
}

inline bool
CXLayer::test_and_add(int64_t x, int64_t y) {
    if (in_use.count(x) || in_use.count(y)) { return false; }
    for (uint64_t q : {x,y}) {
        in_use.insert(q);
        operands.push_back(q);
    }
    return true;
}

inline size_t
CXManager::get_depth() {
    return layers.size();
}

inline void
CXManager::flush(qes::Program<>& program) {
    for (const CXLayer& layer : layers) {
        program.emplace_back("cx", layer.operands);
    }
    layers.clear();
    operand_layer_map.clear();
}

inline qes::Program<>
CXManager::flush() {
    qes::Program<> prog;
    flush(prog);
    return prog;
}


}   // placc
