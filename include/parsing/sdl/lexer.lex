/*
    author: Suhas Vittal
    date:   26 June 2023

    Lexer for "Schedule Description Language" (SDL).
    
    SDL should have the following descriptions:
        (1) Map check to qubit number (makes ASM legal when parsing)
                i.e. decl   X0, 27
                     decl   <check-id>, <qubit>
        (2) Declare micro schedule
                i.e. mus    27
                                h 27
                                cx 27, 1
                                cx 27, 2
                                ...
                     end
        (3) Schedule micro schedules to occur together
                i.e. group  0, 27, 28, 29
                     group  <group>, <check-id>, <check-id>, ...
        (4) Schedule micro schedules to occur after another
                i.e. order  0, 1
                     order  <group-before>, <group-after>
*/

%option noyywrap

%x ASM

%{

#include "parsing/sdl/parser.tab.h"

%}

%%

(?i:decl)           { return DECL; }
(?i:mus)            { BEGIN(ASM); return MUS; }
(?i:end)            { BEGIN(INITIAL); return END; }
[A-Za-z]+           {
                        memcpy(yylval.text, yytext, 128);
                        for (uint i = 0; i < yyleng; i++) {
                            if (yylval.name[i] < 'a') {
                                yylval.name[i] += 'a' - 'A';
                            }
                        }
                        return INST;
                    }
[XZ][0-9]+          {
                        memcpy(yylval.text, yytext, 128);
                        return CHECK;
                    }
[0-9]+              {
                        yylval.id = (uint32_t) atoi(yytext);
                        return ID;
                    }
<ASM>[A-Za-z,0-9]+  {
                        memcpy(yylval.text, yytext, 128);
                        return INST;
                    }
,                   { return SEP; }
[ \t]               { /* ignore whitespace */ }
\n                  { return EOL; }

%%
