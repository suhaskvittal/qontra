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
static std::map<uint32_t, schedule_t>   id_to_schedule;

static std::map<uint32_t, std::set<uint32_t>>
                                        dependents;
static std::map<uint32_t, std::set<uint32_t>>
                                        dependences;

void
sdl_reset_parser() {
    id_to_check.clear();
    check_to_id.clear();
    id_to_schedule.clear();
    dependents.clear();
    dependences.clear();
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
    if (id_to_schedule.count(id)) return -1;
    schedule_t sch = schedule_from_text(std::string(prog.text));
    id_to_schedule[id] = sch;
    if (prog.text != NULL)  free(prog.text);
    return 0;
}

void
sdl_add_dependency(uint32_t id, struct __sdl_ordering ord) {
    for (int i = 0; i < ord.size; i++) {
        dependents[id].insert(ord.dep[i]);
        dependences[ord.dep[i]].insert(id);
    }
    free(ord.dep);
}

namespace qontra {
namespace protean {
namespace compiler {

SDGraph
build_schedule_graph_from_sdl(std::string filename) {
    SDGraph graph;
    sdvertex_t* root = new sdvertex_t;
    root->id = 0;
    graph.add_vertex(root);

    FILE* fin = fopen(filename.c_str(), "r");
    sdl_yystart_file(fin);
    sdl_yyparse_safe();
    fclose(fin);

    std::deque<uint32_t> sat_checks;    // A deque of checks whose dependences
                                        // are already in the graph.
    std::map<uint32_t, std::set<uint32_t>> remaining_dependences;
    for (auto pair : id_to_check) {
        uint32_t ch = pair.second;
        auto dep = dependences[pair.first];
        if (dep.empty())    sat_checks.push_back(ch);
        remaining_dependences[pair.first] = std::set<uint32_t>();
        for (auto x : dep) {
            remaining_dependences[pair.first].insert(x);
        }
    }
    std::map<uint64_t, sdvertex_t*> id_to_vertex;
    while (sat_checks.size()) {
        uint32_t ch = sat_checks.front();
        sat_checks.pop_front();
        
        auto v = new sdvertex_t;
        v->id = (((uint64_t)(ch >> 30)) << ID_TYPE_OFFSET) | (ch & ((1 << 30)-1));

        uint32_t id = check_to_id[ch];
        auto sch = id_to_schedule[check_to_id[ch]];
        // Modify the schedule operands so that any qubit ids
        // are converted to check ids.
        for (auto& inst : sch) {
            for (uint& x : inst.operands) {
                if (x == id) {
                    x = ch;
                }
            }
        }
        v->sch = sch;
        graph.add_vertex(v);
        id_to_vertex[id] = v;
        // Add an edge with the root.
        auto re = new sdedge_t;
        re->src = root;
        re->dst = v;
        re->is_undirected = false;
        graph.add_edge(re);

        for (auto d : dependences[id]) {
            auto w = id_to_vertex[d];
            auto e = new sdedge_t;
            e->src = w;
            e->dst = v;
            e->is_undirected = false;
            graph.add_edge(e);
        }

        for (auto d : dependents[id]) {
            remaining_dependences[d].erase(id);
            if (remaining_dependences[d].empty()) {
                sat_checks.push_back(id_to_check[d]);
            }
        }
    }
    
    return graph;
}

}   // compiler
}   // protean
}   // qontra
