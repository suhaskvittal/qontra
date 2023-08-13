/*
 *  author: Suhas Vittal
 *  date:   5 July 2023
 * */

#ifndef ASM_PARSER_HELPER_h
#define ASM_PARSER_HELPER_h

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// Externs are in parsing/asm.cpp
const extern int    IDLEN;

struct __asm_operand_t {
    uint32_t*   data;
    uint32_t    size;
};

struct __asm_inst_t {
    char                    name[24];
    struct __asm_operand_t  operands;
};

// 512 KB for the program (4K instrutions * 128B).
extern struct __asm_inst_t  ASMParserSchedule[4096];
extern uint32_t             ASMParserScheduleLen;

extern uint64_t     pc;

void        asm_reset_parser(void);

void        asm_clear_labels(void);
void        asm_set_label_pc(char const*, uint64_t);
int         asm_get_label_id(char const*);
int         asm_record_label(char const*);
uint64_t    asm_get_label_pc(int);

#ifdef __cplusplus
}
#endif

#endif  // ASM_PARSER_HELPER_h
