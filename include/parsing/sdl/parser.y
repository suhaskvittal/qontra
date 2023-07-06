/*
    author: Suhas Vittal
    date:   26 June 2023

    Parser for "Schedule Description Language" (SDL).
    See lexer.lex for examples of instructions.
*/

%require "3.2"

%define api.prefix {sdl_yy}

%code requires {

#include "parsing/sdl/common.h"
#include "parsing/sdl/helper.h"

#define __SDL_ERROR(x)  sdl_yyerror(x)

}

%code provides {

int yylex();

extern int  yyparse();
void        yyerror(char const*);

}

%union {
    uint32_t                id;
    char                    check[24];
    struct __sdl_asm_body   prog;
    struct __sdl_ordering   ord;
}

/*
    Tokens:
        DECL = decl (check declaration operation)
        CHECK = a string such as X0 or Z14

        NUM = any number

        MUS = mus (micro-schedule operation)
        ASM = any string: this will be interpreted as ASM
        END = end

        ORD = order (dependency declaration)

        SEP = ,
        EOL = \n
*/

%token  DECL
%token  CHECK 
%token  NUM

%token  MUS
%token  ASM  
%token  END

%token  ORD

%token  SEP 
%token  EOL

%token  ';'
%token  ':'

%type<id>       NUM
%type<check>    CHECK
%type<prog>     ASM
%type<ord>      ordering

%%

program:    /* empty string */
       | DECL CHECK SEP NUM ';' program
{
    if (sdl_declare($4, $2) < 0)    YYABORT;
}
       | MUS NUM ':' ASM END ';' program
{
    if (sdl_assign_schedule($2, $4) < 0)    YYABORT;
}
       | ORD NUM SEP ordering ';' program
{
    sdl_add_dependency($2, $4);
}
       | EOL program
;

ordering:
        NUM SEP ordering
{
    struct __sdl_ordering x;
    x.size = 1 + $3.size;
    x.dep = malloc(x.size * sizeof(uint32_t));
    x.dep[0] = $1;
    memmove(x.dep + 1, $3.dep, $3.size);
    free($3.dep);
    $$ = x;
}
        | NUM
{
    struct __sdl_ordering x;
    x.size = 1;
    x.dep = malloc(1 * sizeof(uint32_t));
    x.dep[0] = $1;
    $$ = x;
}
;

%%

void sdl_yyerror(const char* msg) {
    fprintf(stderr, msg);
}

int sdl_yyparse_safe() {
    return sdl_yyparse();
}
