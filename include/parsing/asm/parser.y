/*
    author: Suhas Vittal
    date:   25 June 2023

    Parser for ASM specified in ISA (see instruction.h)
*/

%require "3.0"

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
    int64_t                 integer;
    double                  decimal;
    char*                   name;
    struct __asm_inst_t     instruction;

    struct __asm_wildcard_data_t anydata;
}

%token ID
%token INTEGER
%token DECIMAL
%token ':'
%token ';'
%token ANNOTATION
%token PROPERTY
%token SEP
%token EOL

%type<integer>      INTEGER
%type<decimal>      DECIMAL
%type<name>         ID

%type<instruction>  instruction
%type<instruction>  operands

%type<anydata> anyval

%%

program: 
        /* empty string */
        | instruction program
{ 
    pc++; 
    asm_add_instruction($1);

    free($1.name);
    if ($1.size > 0) {
        free($1.operands); 
    }
}
        | ANNOTATION ID program
{
    struct __asm_supplement_t annot = { $2 };
    asm_add_annotation(annot);
    free($2);
}
        | PROPERTY ID anyval program
{
    struct __asm_supplement_t prop = { $2 };
    prop.integer = $3.integer;
    prop.decimal = $3.decimal;
    asm_add_property(prop);
}
        | ID ':' program
{
    asm_set_label_inv_pc($1, pc);
    free($1);
}
        | EOL program
;

anyval:
      /* empty */
{
    struct __asm_wildcard_data_t x = {0, 0};
    $$ = x;
}
      | INTEGER anyval
{
    $$ = $2;
    $$.integer = $1;
}
      | DECIMAL anyval
{
    $$ = $2;
    $$.decimal = $1;
}

instruction:
           ID ';'
{
    struct __asm_inst_t inst = DEFAULT_INST;
    inst.name = $1;
    inst.size = 0;
    $$ = inst;
}
           | ID operands ';'
{
    struct __asm_inst_t inst = $2;
    inst.name = $1;
    $$ = inst;
}
;

operands:
        ID
{
    int label = asm_get_label_id($1);
    
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.integer = label;
    x.integer_valid = 1;
    $$ = asm_create_asm_inst_t(x);
}
        | INTEGER
{
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.integer = $1;
    x.integer_valid = 1;
    $$ = asm_create_asm_inst_t(x);
}
        | DECIMAL
{
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.decimal = $1;
    x.decimal_valid = 1;
    $$ = asm_create_asm_inst_t(x);
}
        | ID SEP operands
{
    int label = asm_get_label_id($1);
    
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.integer = label;
    x.integer_valid = 1;

    struct __asm_inst_t inst = $3;
    asm_extend_asm_inst_t(&inst, x);
    $$ = inst;
}
        | INTEGER SEP operands
{
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.integer = $1;
    x.integer_valid = 1;

    struct __asm_inst_t inst = $3;
    asm_extend_asm_inst_t(&inst, x);
    $$ = inst;
}
        | DECIMAL SEP operands
{
    struct __asm_operand_t x = DEFAULT_OPERAND;
    x.decimal = $1;
    x.decimal_valid = 1;

    struct __asm_inst_t inst = $3;
    asm_extend_asm_inst_t(&inst, x);
    $$ = inst;
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

