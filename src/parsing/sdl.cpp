/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "graph/dependence_graph.h"
#include "graph/tanner_graph.h"
#include "instruction.h"
#include "parsing/sdl/common.h"
#include "parsing/sdl/helper.h"

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
    }
    for (auto d : dep)  check_dependences[d].insert(ch);
    free(ord.dep);
}

namespace qontra {
namespace graph {

DependenceGraph
build_dependence_graph_from_sdl(std::string filename) {
    DependenceGraph graph;

    FILE* fin = fopen(filename.c_str(), "r");
    sdl_yystart_file(fin);
    sdl_yyparse_safe();
    fclose(fin);

    std::deque<uint32_t> sat_checks;    // A deque of checks whose dependences
                                        // are already in the graph.
    for (auto pair : check_dependences) {
        if (pair.second.empty())    sat_checks.push_back(pair.first);
    }
    std::map<uint32_t, std::vector<uint32_t>>
        remaining_dependences(check_dependences);
    while (sat_checks.size()) {
        uint32_t ch = sat_checks.front();
        sat_checks.pop_front();
        for (auto& x : check_to_schedule[ch]) {
            Instruction* inst = new Instruction;
            inst->name = x.name;
            inst->operands = x.operands;
            if ((ch >> 30) & tanner::vertex_t::XPARITY) {
                inst->is_measuring_x_check = true;
            }
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

}   // graph
}   // qontra
