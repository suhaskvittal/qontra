/*
    author: Suhas Vittal
    date:   25 June 2023

    Parser for ASM specified in ISA (see instruction.h)
*/

%require "3.2"

%define api.prefix {asm_yy}

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
    char*                   name;
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
    asm_set_label_inv_pc($1, pc);
    free($1);
}
        | EOL program
;

instruction:
           INST ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;
    inst.operands.size = 0;

    asm_add_instruction(inst);

    free(inst.name);
}
           | INST ID ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;

    int label = asm_get_label_id($2);
    inst.operands.size = 1;
    inst.operands.data = malloc(1 * sizeof(uint32_t));
    inst.operands.data[0] = label;

    asm_add_instruction(inst);

    free(inst.name);
    free(inst.operands.data);
}
           | INST ID SEP operands ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;

    int label = asm_get_label_id($2);
    inst.operands.size = 1 + $4.size;
    inst.operands.data = malloc(inst.operands.size * sizeof(uint32_t));
    inst.operands.data[0] = label;
    memmove(inst.operands.data+1, $4.data, $4.size*sizeof(uint32_t));
    free($4.data);

    asm_add_instruction(inst);
    free(inst.name);
    free(inst.operands.data);
}
           | INST operands ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;
    inst.operands.data = $2.data;
    inst.operands.size = $2.size;

    asm_add_instruction(inst);
    free(inst.name);
    free(inst.operands.data);
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
    memcpy(x.data+1, $3.data, $3.size*sizeof(uint32_t));
    free($3.data);
    $$ = x;
}
;

%%

void
asm_yyerror(const char* msg) {
    fprintf(stderr, "asm parsing error: %s\n", msg);
}

int
asm_yyparse_safe() {
    return asm_yyparse();
}

