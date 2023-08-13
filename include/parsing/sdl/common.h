/*
 *  author: Suhas Vittal
 *  date:   26 June 2023
 * */

#ifndef SDL_PARSER_COMMON_h
#define SDL_PARSER_COMMON_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void    sdl_yystart_file(FILE*);
int     sdl_yyparse_safe();    

#ifdef __cplusplus
}
#endif

#endif  // SDL_PARSER_COMMON_h
