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

extern uint64_t pc;

struct __asm_operand_t {
    uint32_t*   data;
    uint32_t    size;
};

struct __asm_inst_t {
    char*                   name;
    struct __asm_operand_t  operands;
};

struct __asm_annotation_t {
    char*   name;
}

void        asm_reset_parser(void);

void        asm_add_instruction(struct __asm_inst_t);
void        asm_add_annotation(struct __asm_annotation_t);

// During parsing, we will replace all label operands with an ID,
// and go back later and replace them with their PC.
//
// A label ID will be negative to indicate it cannot be the PC.

int         asm_declare_label(const char*);
int         asm_set_label_inv_pc(const char*, uint64_t);    
                                            // Inverse pc, as the parser goes
                                            // from bottom to top. We will fix
                                            // this after parsing finishes.
uint64_t    asm_get_label_id(const char*);

#ifdef __cplusplus
}
#endif

#endif  // ASM_PARSER_HELPER_h
