/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "graph/tanner_graph.h"
#include "instruction.h"
#include "parsing/sdl/common.h"
#include "parsing/sdl/helper.h"
#include "protean/compiler.h"

#include <map>

using namespace qontra;

// Parser helper code.

static std::map<uint32_t, uint32_t>     id_to_check;
static std::map<uint32_t, uint32_t>     check_to_id;
static std::map<uint32_t, schedule_t>   check_to_schedule;

static std::map<uint32_t, std::vector<uint32_t>>
                                        check_dependents;

void
reset_parser() {

}

int
declare(uint32_t id, char* check_str) {
    // Convert check string to a integer.
    bool is_x_check = check_str[0] == 'X';
    int n = atoi(check_str+1);
    
    uint32_t ch = n;
    if (is_x_check) ch |= graph::tanner::vertex_t::XPARITY;
    else            ch |= graph::tanner::vertex_t::ZPARITY;

    if (check_to_id.count(ch))  return -1;
    if (id_to_check.count(id))  return -1;
    // So, the ID is what is used in the ASM/SDL code. The
    // check is an integer that we will use everywhere else.
    id_to_check[id] = ch;
    check_to_id[ch] = id;
    return 0;
}

int
assign_schedule(uint32_t id, struct __sdl_asm_body prog) {
    uint32_t check = id_to_check[id];
    schedule_t sch = schedule_from_text(std::string(prog.text));
    check_to_schedule[check] = sch;
}

void
add_dependency(uint32_t id, struct __sdl_ordering ord) {
    std::vector<uint32_t> dep(ord.dep, ord.dep + ord.size);
    check_dependents[id] = dep;
    free(ord.dep);
}

