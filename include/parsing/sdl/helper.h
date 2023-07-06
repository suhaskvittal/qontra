/*
 *  author: Suhas Vittal
 *  date:   5 July 2023
 * */

#ifndef SDL_PARSER_HELPER_h
#define SDL_PARSER_HELPER_h

#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

struct __sdl_asm_body { // Assembly instruction.
    char*       text;
    uint32_t    size;
};

struct __sdl_ordering {
    uint32_t*   dep;
    uint32_t    size;
};

void    reset_parser(void);

int     declare(uint32_t id, char* check_str);
int     assign_schedule(uint32_t check_id, struct __sdl_asm_body);
void    add_dependency(uint32_t, struct __sdl_ordering);

#ifdef __cplusplus
}
#endif

#endif  // SDL_PARSER_HELPER_h
