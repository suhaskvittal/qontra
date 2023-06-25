/*
    author: Suhas Vittal
    date:   25 June 2023

    Parser for ASM specified in ISA (see instruction.h)
*/

%require "3.2"

%code requires {

#include "parsing/asm/common.h"


}

%code provides {

int yylex( );

}

%union {
    uint32_t                arg;
    char                    name[8];
    struct __asm_operand_t  operands;
}

%token INST
%token SEP
%token ARG
%token EOL

%type<arg>      ARG
%type<operands> operands
%type<name>     INST

%%

program: 
        /* empty string */
        | program INST operands EOL
{
    struct __asm_inst_t inst;
    memmove(inst.name, $2, 8);
    memmove(inst.operands.data, $3.data, $3.size);
    inst.operands.size = $3.size;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
;

operands:
        ARG
{ 
    struct __asm_operand_t x;
    x.data[0] = $1;
    x.size = 1;
    $$ = x;
}
        | operands SEP ARG
{ 
    struct __asm_operand_t x = $1;
    x.data[x.size++] = $3;
    $$ = x;
}
;

%%

void
yyerror(const char* msg) {
    fprintf(stderr, "asm parsing error: %s\n", msg);
}

