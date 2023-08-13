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

struct __asm_inst_t     ASMParserSchedule[4096];
uint32_t                ASMParserScheduleLen = 0;

uint64_t    pc = 0;

const int   IDLEN = 24;

// The below data should be hidden from other files.
// For example, if we are using C++ types.

typedef struct {
    int         id;
    uint64_t    pc = 0;
} label_data_t;

static std::map<std::string, label_data_t>    ASMParserLabels; 

// Parser helper code.

void
asm_reset_parser() {
    ASMParserScheduleLen = 0;
    pc = 0;

    asm_clear_labels();
}

void
asm_clear_labels() {
    ASMParserLabels.clear();
}

void
asm_set_label_pc(const char* label, uint64_t x) {
    std::string s(label);
    if (!ASMParserLabels.count(s))  asm_record_label(label);
    ASMParserLabels[s].pc = x;
}

int
asm_get_label_id(const char* label) {
    std::string s(label);
    if (!ASMParserLabels.count(s))  return -1;
    return ASMParserLabels[s].id;
}

int
asm_record_label(const char* label) {
    std::string s(label);
    int id = ASMParserLabels.size();
    ASMParserLabels[s] = {id, 0};
    return id;
}

uint64_t
asm_get_label_pc(int id) {
    for (auto pair : ASMParserLabels) {
        if (pair.second.id == id)   return pair.second.pc;
    }
    return 0;
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
    // Convert C-like schedule to C++.
    schedule_t sch;

    for (uint i = 0; i < ASMParserScheduleLen; i++) {
        __asm_inst_t x = ASMParserSchedule[i];
        std::string name(x.name);
        std::vector<uint> operands;
        for (uint j = 0; j < x.operands.size; j++) {
            operands.push_back(x.operands.data[j]);
        }
        // Check if any operands are labels.
        // 
        // Note that any PCs are inverted, so we must subtract the final
        // pc from the PC in the ASMParserLabels table
        if (ARE_JMP_OR_BR.count(name)) {
            operands[0] = pc - asm_get_label_pc(operands[0]);
        }
        sch.push_back({name, operands, {}});
    }

    return sch;
}

}   // qontra
