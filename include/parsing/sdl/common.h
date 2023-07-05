/*
 *  author: Suhas Vittal
 *  date:   26 June 2023
 * */

#ifndef SDL_PARSER_COMMON_h
#define SDL_PARSER_COMMON_h

#ifdef __cplusplus
extern "C" {
#endif

void    sdl_yystart(FILE*);
int     sdl_yyparse();    

#ifdef __cplusplus
}
#endif

#endif  // SDL_PARSER_COMMON_h
