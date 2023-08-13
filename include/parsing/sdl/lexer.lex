/*
    author: Suhas Vittal
    date:   26 June 2023

    Lexer for "Schedule Description Language" (SDL).
    
    SDL should have the following descriptions:
        (1) Map check to qubit number (makes ASM legal when parsing)
                i.e. decl   X0, 27;
                     decl   <check-id>, <qubit>;
        (2) Declare micro schedule
                i.e. mus    27:
                                h 27;
                                cx 27, 1;
                                cx 27, 2;
                                ...
                     end
        (4) Schedule first mus operand to appear BEFORE the other mus operands
                i.e. order  27, 29;
                     order  <mus-before>, <mus-after1>, <mus-after2>, ...;
*/

%option noyywrap
%option prefix="sdl_yy"

%x ASM

%{

#include "parsing/sdl/parser.tab.h"

#define yylval sdl_yylval

%}

%%

<INITIAL,ASM>\n     { return EOL; }
<ASM>(?i:end)       { BEGIN(INITIAL); return END; }
<ASM>.*;            {
                        yylval.inst = malloc(yyleng+1);
                        memcpy(yylval.inst, yytext, yyleng);
                        yylval.inst[yyleng] = '\0';
                        return INST;
                    }
#.+?\n      { /* this is a comment */ }
(?i:decl)           { return DECL; }
(?i:mus)            { return MUS; }
(?i:order)          { return ORD; }
[XZ][0-9]+          {
                        memcpy(yylval.check, yytext, yyleng);
                        yylval.check[yyleng] = '\0';
                        return CHECK;
                    }
[0-9]+              {
                        yylval.id = atoi(yytext);
                        return NUM;
                    }
,                   { return SEP; }
:                   { BEGIN(ASM); return ':'; }
;                   { return ';'; }
[ \t]               { /* ignore whitespace */ }

%%

void sdl_yystart_file(FILE* fin) {
    sdl_reset_parser();
    sdl_yyrestart(fin);
}
