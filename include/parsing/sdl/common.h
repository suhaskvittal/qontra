/*
 *  author: Suhas Vittal
 *  date:   26 June 2023
 * */

#ifndef SDL_COMMON_h
#define SDL_COMMON_h

#include <stdlib.h>
#include <string.h>

struct __sdl_operand_t {
    uint32_t    data[31];
    uint32_t    size;
};

struct __sdl_body_t {
    char*       inst[128];
    uint32_t    size;
};

struct __sdl_mus_t {
    struct __sdl_body_t sch;
    uint32_t            group = -1;
};

extern char                 SDLParserDeclarations[8][4096];
extern struct __sdl_mus_t   SDLParserSchedules[4096];
extern uint32_t             SDLParserScheduleSize;
extern uint16_t             SDLGroupDependences[16];    // Array of bitvectors

#ifdef __cplusplus
extern "C" {
#endif

void    sdl_yystart(FILE*);
int     sdl_yyparse();    

#ifdef __cplusplus
}
#endif

#endif  // SDL_COMMON_h
