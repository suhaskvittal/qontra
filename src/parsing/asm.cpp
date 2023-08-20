/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "instruction.h"
#include "parsing/asm/common.h"
#include "parsing/asm/helper.h"

#include <map>
#include <string>

// Declarations of external variables in the header.
uint64_t    pc = 0;

// The below data should be hidden from other files.
// For example, if we are using C++ types.

using namespace qontra;

static schedule_t running_schedule;

static std::map<std::string, int> label_to_id;
static std::map<int, uint64_t> label_id_to_pc;

static int label_id_ctr = -1;

// Parser helper code.
void
asm_reset_parser() {
    pc = 0;

    running_schedule.clear();
    label_to_id.clear();
    label_id_to_pc.clear();
}

void
asm_add_instruction(struct __asm_inst_t inst_data) {
    Instruction inst;
    inst.name = std::string(inst_data.name);

    if (std::find(ISA.begin(), ISA.end(), inst.name) == ISA.end()) {
        inst.name = "nop";
    }
    for (uint i = 0; i < inst_data.operands.size; i++) {
        inst.operands.push_back(inst_data.operands.data[i]);
    }
    running_schedule.push_back(inst);
}

void
asm_add_annotation(struct __asm_annotation_t annot) {
    Instruction& inst = running_schedule.back();
    std::string annot_str(annot.name);
    if (annot_str == "no_error") {
        inst.annotations.insert(Annotation::no_error);
    } else if (annot_str == "no_tick") {
        inst.annotations.insert(Annotation::no_tick);
    } else if (annot_str == "round_start") {
        inst.annotations.insert(Annotation::round_start);
    } else if (annot_str == "inject_timing_error") {
        inst.annotations.insert(Annotation::inject_timing_error);
    } else if (annot_str == "flag") {
        inst.annotations.insert(Annotation::flag);
    }
}

int
asm_declare_label(const char* label) {
    std::string s(label);
    if (label_to_id.count(s))   return 1;
    
    label_to_id[s] = label_id_ctr--;
    return 0;
}

int
asm_set_label_inv_pc(const char* label, uint64_t x) {
    std::string s(label);
    if (!label_to_id.count(s))      return 1;
    int id = label_to_id[s];
    if (label_id_to_pc.count(id))   return 1;
    label_id_to_pc[id] = x;
    return 0;
}

uint64_t
asm_get_label_id(const char* label) {
    std::string s(label);
    if (!label_to_id.count(s)) {
        asm_declare_label(label);
    }
    return label_to_id[s];
}

// Glue code from other headers.

namespace qontra {

schedule_t
schedule_from_file(std::string fname) {
    FILE* fin = fopen(fname.c_str(), "r");
    asm_yystart_file(fin);
    asm_yyparse_safe();
    fclose(fin);
    return schedule_after_parse();
}

schedule_t
schedule_from_text(std::string text) {
    asm_yystart_str(text.c_str());
    asm_yyparse_safe();
    return schedule_after_parse();
}

schedule_t
schedule_after_parse() {
    schedule_t sch(running_schedule);
    std::reverse(sch.begin(), sch.end());
    // Convert any label ids in the schedule.
    for (auto& inst : sch) {
        if (IS_BR_TYPE1.count(inst.name)) {
            inst.operands[0] = pc - label_id_to_pc[inst.operands[0]];
        }
    }
    return sch;
}

}   // qontra
