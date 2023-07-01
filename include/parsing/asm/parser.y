/*
    author: Suhas Vittal
    date:   25 June 2023

    Parser for ASM specified in ISA (see instruction.h)
*/

%require "3.2"

%code requires {

#include "parsing/asm/common.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

}

%code provides {

int yylex();

extern int  yyparse();
void        yyerror(char const*);

int64_t     get_label_pc(const char*);

}

%union {
    uint32_t                arg;
    char                    name[12];
    struct __asm_operand_t  operands;
}

%token ID
%token SEP
%token ARG
%token EOL

%type<arg>      ARG
%type<operands> operands
%type<name>     ID

%%

program: 
        /* empty string */
        | instruction program
{ 
    pc++; 
}
        | ID ':' program
{
    struct __asm_label_t x;
    memcpy(x.name, $1, IDLEN);
    x.pc = pc;
    ASMLabelArray[ASMLabelCount++] = x;
}
        | EOL program
;

instruction:
           ID EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);
    inst.operands.size = 0;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | ID ID EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);

    int64_t lpc = get_label_pc($2);
    if (lpc < 0)    YYABORT;
    inst.operands.data[0] = lpc;

    inst.operands.size = 1;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | ID ID operands EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);

    int64_t lpc = get_label_pc($2);
    if (lpc < 0)    YYABORT;
    inst.operands.data[0] = lpc;

    memcpy(inst.operands.data+1, $3.data, $3.size*sizeof(uint32_t));
    inst.operands.size = 1 + $3.size;
    ASMParserSchedule[ASMParserScheduleLen++] = inst;
}
           | ID operands EOL
{
    struct __asm_inst_t inst;
    memcpy(inst.name, $1, IDLEN);
    memcpy(inst.operands.data, $2.data, $2.size*sizeof(uint32_t));
    inst.operands.size = $2.size;
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

int64_t
get_label_pc(const char* name) {
    for (int i = 0; i < ASMLabelCount; i++) {
        if (strcmp(name, ASMLabelArray[i].name) == 0) {
            return ASMLabelArray[i].pc;
        }
    }
    return -1;
}

/*
    Wrapping functions because yy renaming did not work :(
*/

int
asm_yyparse() {
    return yyparse();
}

