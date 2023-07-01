/*
    author: Suhas Vittal
    date:   25 June 2023

    Lexer for ASM specified in ISA (see instruction.h)
*/

%option noyywrap noinput nounput

%{

#include "parsing/asm/parser.tab.h"

#include <stdlib.h>
#include <string.h>

%}

%%

[0-9]+      { 
                yylval.arg = (uint32_t) atoi(yytext); 
                return ARG; 
            }
[A-Za-z]+   { 
                memcpy(yylval.name, yytext, 8); 
                // Force lower case.
                for (int i = 0; i < yyleng; i++) {
                    if (yylval.name[i] < 'a') {
                        yylval.name[i] += 'a' - 'A';
                    }
                }
                return ID; 
            }
[ \t]       { /* ignore whitespace */ }
,           { return SEP; }
:           { return ':'; }
\n          { return EOL; }

%%

void asm_yystart(FILE* fin) {
    ASMParserScheduleLen = 0;
    ASMLabelCount = 0;
    pc = 0;
    yyrestart(fin);
}
