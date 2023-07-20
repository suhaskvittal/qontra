/*
 *  author: Suhas Vittal
 *  date:   19 July 2023
 * */

#include "decoder/restriction.h"

#include <assert.h>

namespace qontra {
namespace decoder {

using namespace graph;

Decoder::result_t
RestrictionDecoder::decode_error(const syndrome_t& syndrome) {
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    timer.clk_start();
    std::vector<uint> detectors = get_nonzero_detectors(syndrome);
    // Get matching results for each restricted lattice.
    auto r_latRG = decode_restricted_lattice(detectors, Color::red, Color::green);
    auto r_latRB = decode_restricted_lattice(detectors, Color::red, Color::blue);
    auto r_latGB = decode_restricted_lattice(detectors, Color::blue, Color::green);
    // Now, we must compute the connected components.
    //
    // First, we just need to figure out which vertices are connected. If two vertices
    // are connected through the boundary, then we separate the connection into
    // two boundary connections.
    std::vector<std::pair<Decoder::result_t*, Color>> results{
        std::make_pair(&r_latRG, Color::blue),
        std::make_pair(&r_latRB, Color::green),
        std::make_pair(&r_latGB, Color::red)
    };
    typedef std::pair<decoding::vertex_t*, Color>   colored_vertex_t;
    std::map<colored_vertex_t, std::set<colored_vertex_t>> connections;

    for (auto pair : results) {
        Decoder::result* res = pair.first;
        Color restricted_color = pair.second;
        for (auto aa : res->error_assignments) {
            auto v1 = decoding_graph.get_vertex(d1);
            auto v2 = decoding_graph.get_vertex(d2);
            if (v1->id == BOUNDARY_INDEX)   std::swap(v1, v2);
            auto path = decoding_graph.get_error_chain_data(v1, v2).error_chain;

            std::vector<Color> colors_along_path(path.size());
            colors_along_path[0] = detector_to_color[v1];

            bool found_boundary = false;
            
            colored_vertex_t cv1 = std::make_pair(v1, detector_to_color[v1]);
            colored_vertex_t cv2 = std::make_pair(v1, detector_to_color[v2]);

            for (uint i = 1; i < path.size(); i++) {
                auto vx = path[i-1];
                auto vy = path[i];

                Color yc = Color::red;
                if ((colors_along_path[i-1] == Color::red && restricted_color == Color::blue)
                    || (colors_along_path[i-1] == Color::blue && restricted_color == Color::red))
                {
                    yc = Color::green;
                } else if ((colors_along_path[i-1] == Color::red && restricted_color == Color::green)
                    || (colors_along_path[i-1] == Color::green && restricted_color == Color::red))
                {
                    yc = Color::blue;
                } else if ((colors_along_path[i-1] == Color::green && restricted_color == Color::blue)
                    || (colors_along_path[i-1] == Color::blue && restricted_color == Color::green))
                {
                    yc = Color::red;
                }

                // Assertion that color is correct:
                if (vy->id != BOUNDARY_INDEX) {
                    std::string alert = std::string("colors do not match for D") 
                                        + std::to_string(vy->id) << "!\n";
                    assert(detector_to_color[vy] == yc, alert.c_str());
                } else {
                    found_boundary = true;
                    colored_vertex_t cb = std::make_pair(vy, yc);
                    // Add connection to v1.
                    connections[cb].insert(cv1);
                    connections[cv1].insert(cb);
                    // If v2 is not the boundary, add a connection there as well.
                    if (v2 != vy) {
                        connections[cb].insert(cv2);
                        connections[cv2].insert(cb);
                    }
                }
                colors_along_path[i] = yc;
            }
            // If the boundary was not found, then just connect v1 to v2. 
            if (!found_boundary) {
                connections[cv1].insert(cv2);
                connections[cv2].insert(cv1);
            }
        }
    }
    // Compute all connected components originating from the red boundary.
    //
    // We will just do a DFS for this.
}

}   // decoder
}   // qontra
