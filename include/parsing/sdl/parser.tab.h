/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_SDL_PARSER_TAB_H_INCLUDED
# define YY_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_SDL_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif
/* "%code requires" blocks.  */
#line 11 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.y"


#include "parsing/sdl/common.h"
#include "parsing/sdl/helper.h"

#define __SDL_ERROR(x)  sdl_yyerror(x)


#line 58 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.tab.h"

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    DECL = 258,                    /* DECL  */
    CHECK = 259,                   /* CHECK  */
    NUM = 260,                     /* NUM  */
    MUS = 261,                     /* MUS  */
    ASM = 262,                     /* ASM  */
    END = 263,                     /* END  */
    ORD = 264,                     /* ORD  */
    SEP = 265,                     /* SEP  */
    EOL = 266                      /* EOL  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 29 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.y"

    uint32_t                id;
    char                    check[24];
    struct __sdl_asm_body   prog;
    struct __sdl_ordering   ord;

#line 93 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);

/* "%code provides" blocks.  */
#line 20 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.y"


int yylex();

extern int  yyparse();
void        yyerror(char const*);


#line 117 "/Users/svittal/Documents/qontra/include/parsing/sdl/parser.tab.h"

#endif /* !YY_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_SDL_PARSER_TAB_H_INCLUDED  */
