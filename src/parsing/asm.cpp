/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "instruction.h"
#include "parsing/asm/common.h"

#include <map>
#include <string>

// Declarations of external variables in the header.

struct __asm_inst_t     ASMParserSchedule[4096];
uint32_t                ASMParserScheduleLen = 0;

uint64_t    pc = 0;

const int   IDLEN = 24;
const int   MAX_OPERANDS = 25;

// The below data should be hidden from other files.
// For example, if we are using C++ types.

typedef struct {
    int         id;
    uint64_t    pc = 0;
} label_data_t;

static std::map<std::string, label_data_t>    ASMParserLabels; 

// Parser helper code.

void
reset_parser() {
    ASMParserScheduleLen = 0;
    pc = 0;

    clear_labels();
}

void
clear_labels() {
    ASMParserLabels.clear();
}

void
set_label_pc(const char* label, uint64_t x) {
    std::string s(label);
    if (!ASMParserLabels.count(s))  record_label(label);
    ASMParserLabels[s].pc = x;
}

int
get_label_id(const char* label) {
    std::string s(label);
    if (!ASMParserLabels.count(s))  return -1;
    return ASMParserLabels[s].id;
}

int
record_label(const char* label) {
    std::string s(label);
    int id = ASMParserLabels.size();
    ASMParserLabels[s] = {id, 0};
    return id;
}

uint64_t
get_label_pc(int id) {
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
    asm_yystart(fin);
    asm_yyparse();
    fclose(fin);
    // Convert C-like schedule to C++.
    schedule_t sch;

    for (uint i = 0; i < ASMParserScheduleLen; i++) {
        __asm_inst_t x = ASMParserSchedule[i];
        std::string name(x.name);
        std::vector<uint> operands;
        operands.assign(x.operands.data, x.operands.data + x.operands.size);
        // Check if any operands are labels.
        // 
        // Note that any PCs are inverted, so we must subtract the final
        // pc from the PC in the ASMParserLabels table
        if (ARE_JMP_OR_BR.count(name)) {
            operands[0] = pc - get_label_pc(operands[0]);
        }
        sch.push_back({name, operands, {}});
    }

    return sch;
}

}   // qontra
