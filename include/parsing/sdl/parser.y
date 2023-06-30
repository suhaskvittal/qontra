/*
    author: Suhas Vittal
    date:   26 June 2023

    Parser for "Schedule Description Language" (SDL).
    See lexer.lex for examples of instructions.
*/

%require "3.2"

%code requires {

#include "parsing/sdl/common.h"

#define __SDL_ERROR(x)  sdl_yyerror(x)
                                                        // [x][y] = 1 if x precedes y
}

%code provides {

int yylex();

extern int  yyparse();
void        yyerror(char const*);

}

%union {
    uint32_t                id;
    char                    text[128];
    struct __sdl_operand_t  operands;
    struct __sdl_body_t     body;
}

%token DECL
%token MUS
%token INST
%token CHECK
%token ID
%token SEP
%token EOL
%token END

%type<id>       ID
%type<text>     CHECK
%type<text>     INST
%type<body>     body
%type<operands> operands

%%

program:    /* empty string */
       | instruction EOL program
;

instruction:
           DECL CHECK SEP ID
{
    int len = strlen($2);
    memcpy(SDLParserDeclaration[$4], $2, len+1);
}
           | MUS ID EOL body
{
    for (int i = 0; i < $4.size; i++) {
        SDLParserSchedules[$2].sch.inst[i] = $4.inst[i];
    } 
    SDLParserSchedules[$2].sch.size = $4.size;
}
           | INST ID SEP operands
{
    if (strcmp($1, "group")) {
        uint32_t grp = $2;
        for (int i = 0; i < $4.size; i++) {
            uint32_t id = $4.data[i];
            SDLParserSchedules[id].group = grp;
        }
    } else if (strcmp($1, "order")) {
        uint32_t x = $2;
        uint32_t y = $4.data[0];
        SDLGroupDependences[x] |= (1 << y);
    }
}
;

body:
    INST EOL body
{
    struct __sdl_body_t x;
    int len = strlen($1);
    x.inst[0] = malloc(sizeof(char) * (len+1));
    memcpy(x.inst[0], $1, len+1);
    // Copy body ($3) to x.
    memcpy(x.inst + 1, $3.inst, sizeof(char*) * $3.size);
    x.size = 1 + $3.size;
    $$ = x;
}
    | END
{
    struct __sdl_body_t x;
    x.size = 0;
    $$ = x;
}
;

operands:
        ID SEP operands
{
    struct __sdl_operand_t x;
    x.data[0] = $1;
    memcpy(x.data + 1, $3.data, $3.size);
    x.size = 1 + $3.size;
    $$ = x;
}
        | ID
{
    struct __sdl_operand_t x;
    x.data[0] = $1;
    x.size = 1;
    $$ = x;
}
;

%%

void yyparse(const char* msg) {
    fprintf(stderr, msg);
}

int sdl_yyparse() {
    return yyparse();
}
