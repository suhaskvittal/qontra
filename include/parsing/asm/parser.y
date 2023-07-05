/*
    author: Suhas Vittal
    date:   25 June 2023

    Parser for ASM specified in ISA (see instruction.h)
*/

%require "3.2"

%code requires {

#include "parsing/asm/common.h"
#include "parsing/asm/helper.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

}

%code provides {

int yylex();

extern int  yyparse();
void        yyerror(char const*);

}

%union {
    uint32_t                arg;
    char                    name[24];
    struct __asm_operand_t  operands;
}

/*
    Tokens:
        INST = instructions
        ID   = text that is not an instruction (i.e. LABEL)
        NUM  = numbers
        SEP  = separator (,)
*/

%token INST
%token ID
%token NUM
%token ':'
%token ';'
%token SEP
%token EOL

%type<arg>          NUM
%type<name>         INST
%type<name>         ID
%type<operands>     operands

%%

program: 
        /* empty string */
        | instruction program
{ 
    pc++; 
}
        | ID ':' program
{
    set_label_pc($1, pc);
}
        | EOL program
;

instruction:
           INST ';'
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);
    inst.operands.size = 0;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | INST ID ';'
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);

    int label = get_label_id($2);
    if (label < 0) {
        label = record_label($2);    
    }
    inst.operands.size = 1;
    inst.operands.data = malloc(1 * sizeof(uint32_t));
    inst.operands.data[0] = label;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | INST ID operands ';'
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);

    int label = get_label_id($2);
    if (label < 0) {
        label = record_label($2);    
    }
    inst.operands.size = 1 + $3.size;
    inst.operands.data = malloc(inst.operands.size * sizeof(uint32_t));
    inst.operands.data[0] = label;
    memmove(inst.operands.data+1, $3.data, $3.size*sizeof(uint32_t));
    free($3.data);
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | INST operands ';'
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);
    inst.operands.data = $2.data;
    inst.operands.size = $2.size;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
;

operands:
        NUM
{
    struct __asm_operand_t x;
    x.data = malloc(1 * sizeof(uint32_t));
    x.data[0] = $1;
    x.size = 1;
    $$ = x;
}
        | NUM SEP operands
{
    struct __asm_operand_t x;
    x.size = 1 + $3.size;
    x.data = malloc(x.size * sizeof(uint32_t));
    x.data[0] = $1;
    memmove(x.data+1, $3.data, $3.size*sizeof(uint32_t));
    free($3.data);
    $$ = x;
}
;

%%

void
yyerror(const char* msg) {
    fprintf(stderr, "asm parsing error: %s\n", msg);
}

/*
    Wrapping functions because yy renaming did not work :(
*/

int
asm_yyparse() {
    int parse_out = yyparse();
    // Free any heap-allocated memory.
    for (int i = 0; i < ASMParserScheduleLen; i++) {
        struct __asm_inst_t* inst = &ASMParserSchedule[i];
        if (inst->operands.size) free(inst->operands.data);
    }
    return parse_out;
}

