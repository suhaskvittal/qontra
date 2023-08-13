/*
    author: Suhas Vittal
    date:   25 June 2023

    Lexer for ASM specified in ISA (see instruction.h)
*/

%option noyywrap noinput nounput

%option prefix="asm_yy"

%{

#include "parsing/asm/parser.tab.h"

#include <stdlib.h>
#include <string.h>

void    force_lowercase(char*, int);

#define yylval  asm_yylval

%}

%%

#.+?\n      { /* this is a comment */ }

[0-9]+      { 
                yylval.arg = (uint32_t) atoi(yytext); 
                return NUM; 
            }
[A-Za-z]+   {
                memcpy(yylval.name, yytext, 24);
                force_lowercase(yylval.name, yyleng);
                return INST;
            }
[A-Za-z_0-9]+  { 
                memcpy(yylval.name, yytext, 24); 
                return ID; 
            }
[ \t]       { /* ignore whitespace */ }
:           { return ':'; }
;           { return ';'; }
,           { return SEP; }
\n          { BEGIN(INITIAL); return EOL; }

%%

void force_lowercase(char* text, int len) {
    for (int i = 0; i < len; i++) {
        if (text[i] < 'a') {
            text[i] += 'a' - 'A';
        }
    }
}

void asm_yystart_file(FILE* fin) {
    asm_reset_parser();
    asm_yyrestart(fin);
}

void asm_yystart_str(const char* text) {
    asm_reset_parser();
    asm_yy_scan_string(text);
}

