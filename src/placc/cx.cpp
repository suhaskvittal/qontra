/*
 *  author: Suhas Vittal
 *  date:   30 May 2024
 * */

#include "placc/cx.h"

namespace placc {

CXManager::CXManager()
    :layers(),
    operand_layer_map()
{}

void
CXManager::push_back_cx(uint64_t x, uint64_t y) {
    size_t xstart = operand_layer_map.count(x) ? operand_layer_map.at(x)+1 : 0;
    size_t ystart = operand_layer_map.count(y) ? operand_layer_map.at(y)+1 : 0;
    for (size_t i = std::max(xstart,ystart); i < layers.size(); i++) {
        if (layers[i].test_and_add(x,y)) {
            operand_layer_map[x] = i;
            operand_layer_map[y] = i;
            return;
        }
    }
    layers.emplace_back();
    layers.back().test_and_add(x,y);
    operand_layer_map[x] = layers.size()-1;
    operand_layer_map[y] = layers.size()-1;
}

}   // placc
