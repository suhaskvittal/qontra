/*
    author: Suhas Vittal
    date:   22 June 2023

    Simple schedule lexer using a restricted version of QASM.
*/

%option noyywrap
%option c++
%option yyclass="ProteanScheduleLexer"

%x      QREG
%x      GATE

%{

#include "common.h"

#include "v1.tab.hh"

char    qreg[1024];
int     qreg_len;
int     n_qubits;

%}

%%

qreg                            { BEGIN(QREG); return QREG; }
<QREG,GATE>[A-Za-z][A-Za-z0-9]+ { yylval = yytext; return NAME; }
[A-Za-z][A-Za-z0-9]+            {   
                                    yylval = yytext; 
                                    for (uint i = 0; i < strlen(yytext); i++) {
                                        // Force uppercase.
                                        if (yylval[i] > 'Z')    yylval[i] -= 'A';
                                    }
                                    BEGIN(GATE); 
                                    return GATE; 
                                }
\d+                             { yylval = atoi(yytext); return ARG; }
\[                              { return '['; }
\]                              { return ']'; }
,                               { return SEP; }
;                               { return EOL; }
$done                           { return DONE; }

%%
