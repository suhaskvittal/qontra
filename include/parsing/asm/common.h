/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#ifndef ASM_PARSER_LEXER_COMMON_h
#define ASM_PARSER_LEXER_COMMON_h

#include <stdint.h>
#include <stdio.h>
#include <string.h>

extern "C"  void asm_yystart(FILE*);
extern "C"  int yyparse();
extern "C"  void yyerror(char const*);
extern "C"  int yylex();

#include "parsing/asm/parser.tab.h"
#include "parsing/asm/lex.tab.h"

#endif  // ASM_PARSER_COMMON_h
