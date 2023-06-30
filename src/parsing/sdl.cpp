/*
 *  author: Suhas Vittal
 *  date:   25 June 2023
 * */

#include "parsing/sdl/common.h"

char                SDLParserDeclarations[8][4096];
struct __sdl_mus_t  SDLParserSchedules[4096];
uint32_t            SDLParserScheduleSize = 0;
uint16_t            SDLGroupDependences[16];

