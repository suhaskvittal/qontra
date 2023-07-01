/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#ifndef ASM_PARSER_LEXER_COMMON_h
#define ASM_PARSER_LEXER_COMMON_h

#include <stdint.h>
#include <stdio.h>
#include <string.h>

struct __asm_operand_t {
    uint32_t    data[28];
    uint32_t    size;
};

struct __asm_inst_t {   // Each instruction is 128B.
    char                    name[12];   // 12 B
    struct __asm_operand_t  operands;   // 116 B
};

// 512 KB for the program (4K instrutions * 128B).
extern struct __asm_inst_t  ASMParserSchedule[4096];
extern uint32_t             ASMParserScheduleLen;

#ifdef __cplusplus
extern "C" {
#endif

void asm_yystart(FILE*);
int  asm_yyparse();

#ifdef __cplusplus
}
#endif

#endif  // ASM_PARSER_COMMON_h
