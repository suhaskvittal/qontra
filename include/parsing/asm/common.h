/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#ifndef ASM_COMMON_h
#define ASM_COMMON_h

#include <stdint.h>
#include <stdio.h>
#include <string.h>

struct __asm_operand_t {
    uint32_t    data[31];
    uint32_t    size;
};

struct __asm_inst_t {   // Each instruction is 128B.
    char                    name[8];    // 8 B
    struct __asm_operand_t  operands;   // 120 B
};

// 512 KB for the program (4K instrutions * 128B).
extern struct __asm_inst_t  ASMParserSchedule[4096];
extern uint32_t             ASMParserScheduleLen;

#ifdef  __cplusplus
extern "C" {
#endif

/*
void        asm_yystart(FILE*);
extern int  asm_yyparse();
void        asm_yyerror(char const*);
*/
void        asm_yystart(FILE*);
int         asm_yyparse();
void        asm_yyerror(char const*);

#ifdef  __cplusplus
}   // extern "C"
#endif

#endif  // ASM_PARSER_COMMON_h
