/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "graph/tanner_graph.h"
#include "instruction.h"
#include "parsing/sdl/common.h"
#include "parsing/sdl/helper.h"
#include "protean/compiler.h"

#include <deque>
#include <map>
#include <set>
#include <vector>

using namespace qontra;

// Parser helper code.

static std::map<uint32_t, uint32_t>     id_to_check;
static std::map<uint32_t, uint32_t>     check_to_id;
static std::map<uint32_t, schedule_t>   check_to_schedule;

static std::map<uint32_t, std::set<uint32_t>>
                                        check_dependents;
static std::map<uint32_t, std::set<uint32_t>>
                                        check_dependences;

void
sdl_reset_parser() {

}

int
sdl_declare(uint32_t id, char* check_str) {
    // Convert check string to a integer.
    bool is_x_check = check_str[0] == 'X';
    int n = atoi(check_str+1);
    
    uint32_t ch = n;
    if (is_x_check) ch |= graph::tanner::vertex_t::XPARITY << 30;
    else            ch |= graph::tanner::vertex_t::ZPARITY << 30;

    if (check_to_id.count(ch))  return -1;
    if (id_to_check.count(id))  return -1;
    // So, the ID is what is used in the ASM/SDL code. The
    // check is an integer that we will use everywhere else.
    id_to_check[id] = ch;
    check_to_id[ch] = id;
    return 0;
}

int
sdl_assign_schedule(uint32_t id, struct __sdl_asm_body prog) {
    uint32_t check = id_to_check[id];
    if (check_to_schedule.count(check)) return -1;
    schedule_t sch = schedule_from_text(std::string(prog.text));
    check_to_schedule[check] = sch;
    free(prog.text);
    return 0;
}

void
sdl_add_dependency(uint32_t id, struct __sdl_ordering ord) {
    uint32_t ch = id_to_check[id];
    for (int i = 0; i < ord.size; i++) {
        check_dependents[ch].insert(id_to_check[ord.dep[i]]);
        check_dependences[id_to_check[ord.dep[i]]].insert(ch);
    }
    free(ord.dep);
}

namespace qontra {
namespace protean {
namespace compiler {

SDGraph
build_dependence_graph_from_sdl(std::string filename) {
    SDGraph graph;

    FILE* fin = fopen(filename.c_str(), "r");
    sdl_yystart_file(fin);
    sdl_yyparse_safe();
    fclose(fin);

    std::deque<uint32_t> sat_checks;    // A deque of checks whose dependences
                                        // are already in the graph.
    for (auto pair : check_dependences) {
        if (pair.second.empty())    sat_checks.push_back(pair.first);
    }
    std::map<uint32_t, std::set<uint32_t>>
        remaining_dependences(check_dependences);
    uint64_t vid = 0;
    while (sat_checks.size()) {
        uint32_t ch = sat_checks.front();
        sat_checks.pop_front();
        
        auto v = new sdvertex_t;
        v->id = ch;
        auto sch = check_to_schedule[ch];
        // Modify the schedule operands so that any qubit ids
        // are converted to check ids.
        for (auto& inst : sch) {
            for (uint& x : inst.operands) {
                if (x == check_to_id[ch]) {
                    x = ch;
                }
            }
        }
        graph.add_vertex(v);

        for (auto d : check_dependences[ch]) {
            auto w = graph.get_vertex(d);
            auto e = new sdedge_t;
            e->src = w;
            e->dst = v;
            e->is_undirected = false;
            graph.add_edge(e);
        }

        for (auto d : check_dependents[ch]) {
            remaining_dependences[d].erase(ch);            
            if (remaining_dependences[d].empty()) {
                sat_checks.push_back(d);
            }
        }
    }
    
    return graph;
}

}   // compiler
}   // protean
}   // qontra
