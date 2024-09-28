/*
 *  author: Suhas Vittal
 *  date:   23 May 2024
 * */

#ifndef CC_GEN_h
#define CC_GEN_h

#include <qontra/graph/tanner_graph.h>
#include <qontra/graph/io.h>
#include <qontra/ext/stim.h>

#include <vtils/filesystem.h>

#include <unordered_map>

// It's just a header file for commonly used main functions. I don't care
// about convention here.
using namespace qontra;
using namespace graph;
using namespace vtils;

DetailedStimCircuit make_capacity(
        TannerGraph*, 
        fp_t,
        bool is_mx,
        const std::unordered_map<sptr<tanner::vertex_t>, int>& color_map);

#include "gen.inl"

#endif  // CC_GEN_h
