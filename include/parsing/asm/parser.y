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

void    yyerror(char const*);
int     yyparse( );

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
        | instruction program
;

instruction:
           INST EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, 8);
    inst.operands.size = 0;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | INST ARG EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, 8);
    inst.operands.data[0] = $2;
    inst.operands.size = 1;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | INST ARG SEP operands EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, 8);
    inst.operands.data[0] = $2;
    memcpy(inst.operands.data+1, $4.data, $4.size*sizeof(uint32_t));
    inst.operands.size = 1 + $4.size;
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
        | ARG SEP operands
{
    struct __asm_operand_t x;
    x.data[0] = $1;
    memcpy(x.data+1, $3.data, $3.size*sizeof(uint32_t));
    x.size = 1 + $3.size;
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
void
asm_yyerror(const char* msg) {
    return yyerror(msg);
}

int
asm_yyparse() {
    return yyparse();
}

