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

struct __asm_operand_t DEFAULT_OPERAND = {0, 0, 0, 0};
struct __asm_inst_t DEFAULT_INST = {NULL, NULL, 0};

// The below data should be hidden from other files.
// For example, if we are using C++ types.

using namespace qontra;

static schedule_t running_schedule;

static std::map<std::string, int> label_to_id;
static std::map<int, uint64_t> label_id_to_pc;

static int label_id_ctr = -1;

// Parser helper code.
struct __asm_inst_t
asm_create_asm_inst_t(struct __asm_operand_t op) {
    struct __asm_inst_t inst;
    inst.operands = (struct __asm_operand_t*)malloc(1 * sizeof(struct __asm_operand_t));
    inst.operands[0] = op;
    inst.size = 1;
    return inst;
}

void
asm_extend_asm_inst_t(struct __asm_inst_t* inst_p, struct __asm_operand_t op) {
    if (inst_p->size == 0) {
        *inst_p = asm_create_asm_inst_t(op);
    } else {
        inst_p->operands = (struct __asm_operand_t*)
                                realloc(inst_p->operands, (inst_p->size+1) * sizeof(struct __asm_operand_t));
        memmove(inst_p->operands+1, inst_p->operands, inst_p->size*sizeof(struct __asm_operand_t));
        inst_p->operands[0] = op;
        inst_p->size++;
    }
}

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
        std::cout << "(asm) could not find instruction " << inst.name << "\n";
        inst.name = "nop";
    }
    if (IS_QUANTUM.count(inst.name) || IS_QUANTUM_LIKE.count(inst.name)) {
        // We can assume all operands are integers.
        for (uint i = 0; i < inst_data.size; i++) {
            inst.operands.qubits.push_back(inst_data.operands[i].integer);
        }
    } else if (IS_CLASSICALLY_CONTROLLED.count(inst.name)) {
        for (uint i = 0; i < inst_data.size; i += 2) {
            inst.operands.measurements.push_back(inst_data.operands[i].integer);
            inst.operands.qubits.push_back(inst_data.operands[i+1].integer);
        }
    } else if (INSTRUCTION_USES_ANGLES.count(inst.name)) {
        // Now, we need to check for floats, as those are angles.
        for (uint i = 0; i < inst_data.size; i++) {
            if (inst_data.operands[i].integer_valid) {
                inst.operands.qubits.push_back(inst_data.operands[i].integer);
            } else if (inst_data.operands[i].decimal_valid) {
                inst.operands.angles.push_back(inst_data.operands[i].decimal);
            }
        }
    } else if (IS_DECODE_TYPE1.count(inst.name)) {
        inst.operands.events.push_back(inst_data.operands[0].integer);
        for (uint i = 1; i < inst_data.size; i++) {
            inst.operands.measurements.push_back(inst_data.operands[i].integer);
        }
    } else if (IS_DECODE_TYPE2.count(inst.name)) {
        inst.operands.observables.push_back(inst_data.operands[0].integer);
        inst.operands.events.push_back(inst_data.operands[0].integer);
        for (uint i = 1; i < inst_data.size; i++) {
            inst.operands.measurements.push_back(inst_data.operands[i].integer);
        }
    } else if (IS_DECODE_TYPE3.count(inst.name)) {
        inst.operands.observables.push_back(inst_data.operands[0].integer);
        inst.operands.frames.push_back(inst_data.operands[1].integer);
    } else if (IS_DECODE_TYPE4.count(inst.name)) {
        inst.operands.frames.push_back(inst_data.operands[0].integer);
    } else if (IS_BR_TYPE1.count(inst.name)) {
        inst.operands.labels.push_back(inst_data.operands[0].integer);
        inst.operands.events.push_back(inst_data.operands[1].integer);
    }
    running_schedule.push_back(inst);
}

void
asm_add_annotation(struct __asm_supplement_t annot) {
    Instruction& inst = running_schedule.back();
    std::string annot_str(annot.name);
    inst.annotations.insert(annot_str);
}

void
asm_add_property(struct __asm_supplement_t prop) {
    Instruction& inst = running_schedule.back();
    std::string prop_str(prop.name);
    inst.properties[prop_str] = {
        prop.integer,
        prop.decimal
    };
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
        if (IS_BR.count(inst.name)) {
            for (uint i = 0; i < inst.operands.labels.size(); i++) {
                inst.operands.labels[i] = pc - label_id_to_pc[inst.operands.labels[i]];
            }
        }
    }
    return sch;
}

}   // qontra
