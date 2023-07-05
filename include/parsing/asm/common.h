/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#ifndef ASM_PARSER_LEXER_COMMON_h
#define ASM_PARSER_LEXER_COMMON_h

#ifdef __cplusplus
extern "C" {
#endif

void asm_yystart(FILE*);
int  asm_yyparse();

#ifdef __cplusplus
}
#endif

#endif  // ASM_PARSER_COMMON_h
