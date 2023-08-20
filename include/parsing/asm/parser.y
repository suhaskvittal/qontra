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
    struct __asm_inst_t     instruction;
}

/*
    Tokens:
        ID   = text
        NUM  = numbers
        SEP  = separator (,)
*/

%token ID
%token NUM
%token ':'
%token ';'
%token ANNOTATION
%token SEP
%token EOL

%type<arg>          NUM
%type<name>         ID

%type<instruction>  instruction
%type<operands>     operands

%%

program: 
        /* empty string */
        | instruction program
{ 
    pc++; 
    asm_add_instruction($1);

    free($1.name);
    if ($1.operands.size > 0) {
        free($1.operands.data); 
    }
}
        | ANNOTATION ID program
{
    struct __asm_annotation_t annot = { $2 };
    asm_add_annotation(annot);
    free($2);
}
        | ID ':' program
{
    asm_set_label_inv_pc($1, pc);
    free($1);
}
        | EOL program
;

instruction:
           ID ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;
    inst.operands.size = 0;

    $$ = inst;
}
           | ID ID ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;

    int label = asm_get_label_id($2);
    inst.operands.size = 1;
    inst.operands.data = malloc(1 * sizeof(uint32_t));
    inst.operands.data[0] = label;

    $$ = inst;
}
           | ID ID SEP operands ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;

    int label = asm_get_label_id($2);
    inst.operands.size = 1 + $4.size;
    inst.operands.data = malloc(inst.operands.size * sizeof(uint32_t));
    inst.operands.data[0] = label;
    memmove(inst.operands.data+1, $4.data, $4.size*sizeof(uint32_t));
    free($4.data);

    $$ = inst;
}
           | ID operands ';'
{
    struct __asm_inst_t inst;
    inst.name = $1;
    inst.operands.data = $2.data;
    inst.operands.size = $2.size;

    $$ = inst;
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

