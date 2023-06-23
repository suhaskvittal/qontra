/* author: Suhas Vittal
    date:   22 June 2023

    Simple schedule parser using a restricted version of QASM.
*/

%require "3.2"
%language "c++"
%define api.parser.class {ProteanScheduleParser}

%{

#include "common.h"

#include "instruction.h"

using namespace qontra;

schedule_t<qc::Instruction> schedule;

%}

// Tokens:
%token QREG
%token NAME
%token GATE
%token SEP
%token EOL
%token ARG
%token '['
%token ']'

%%

program: /* empty string */
    | decl EOL program
    | op EOL program

decl: 
    QREG NAME '[' ARG ']' EOL
{
    qreg_len = strlen($2);
    memcpy(qreg, $2, qreg_len+1); 
    n_qubits = $4;
};

exp:
   /* empty string */
{ $$ = std::vector<uint>(); }
   | exp SEP NAME '[' ARG ']'
{ if (strcmp($3, qreg) == 0)    $$.push_back($5); }

op: 
  GATE exp
{
    qc::Instruction x;
    x.name = $1;
    x.operands = $2;
    schedule.push_back(x);
}

%%

yyerror(char* m) {
    std::cerr << "parser error: " << m << "\n";
}
