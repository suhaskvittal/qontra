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

%x ASM

%{

#include "parsing/sdl/parser.tab.h"

%}

%%

#.+?\n      { /* this is a comment */ }
(?i:decl)           { return DECL; }
(?i:mus)            { BEGIN(ASM); return MUS; }
(?i:end)            { BEGIN(INITIAL); return END; }
(?i:order)          { return ORD; }
[XZ][0-9]+          {
                        memcpy(yylval.check, yytext, yyleng);
                        return CHECK;
                    }
[0-9]+              {
                        yylval.id = atoi(yytext);
                        return NUM;
                    }
<ASM>.*?            {
                        yylval.prog.text = malloc(yyleng+1);
                        memcpy(yylval.prog.text, yytext, yyleng);
                        yylval.prog.text[yyleng] = '\0';
                        return ASM;
                    }
,                   { return SEP; }
:                   { return ':'; }
;                   { return ';'; }
[ \t]               { /* ignore whitespace */ }
\n                  { return EOL; }

%%

void sdl_yystart(FILE* fin) {
    yyrestart(fin);
}
