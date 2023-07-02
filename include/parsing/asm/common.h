/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#ifndef ASM_PARSER_LEXER_COMMON_h
#define ASM_PARSER_LEXER_COMMON_h

#include <stdint.h>
#include <stdio.h>
#include <string.h>

// Externs are in parsing/asm.cpp
const extern int    IDLEN;
const extern int    MAX_OPERANDS;

struct __asm_operand_t {
    uint32_t    data[25];
    uint32_t    size;
};

struct __asm_inst_t {   // Each instruction is 128B.
    char                    name[24];   // 24 B
    struct __asm_operand_t  operands;   // 104 B
};

// 512 KB for the program (4K instrutions * 128B).
extern struct __asm_inst_t  ASMParserSchedule[4096];
extern uint32_t             ASMParserScheduleLen;

extern uint64_t     pc;

#ifdef __cplusplus
extern "C" {
#endif

void asm_yystart(FILE*);
int  asm_yyparse();

void        reset_parser(void);

void        clear_labels(void);
void        set_label_pc(char const*, uint64_t);
int         get_label_id(char const*);
int         record_label(char const*);
uint64_t    get_label_pc(int);

#ifdef __cplusplus
}
#endif

#endif  // ASM_PARSER_COMMON_h
