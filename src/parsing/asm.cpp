/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "instruction.h"
#include "parsing/asm/common.h"

struct __asm_inst_t     ASMParserSchedule[4096];
uint32_t                ASMParserScheduleLen = 0;

namespace qontra {

schedule_t
from_file(std::string fname) {
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
        sch.push_back({name, operands, {}});
    }

    return sch;
}

}   // qontra
