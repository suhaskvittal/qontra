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

#ifndef YY_ASM_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_ASM_PARSER_TAB_H_INCLUDED
# define YY_ASM_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_ASM_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef ASM_YYDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define ASM_YYDEBUG 1
#  else
#   define ASM_YYDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define ASM_YYDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined ASM_YYDEBUG */
#if ASM_YYDEBUG
extern int asm_yydebug;
#endif
/* "%code requires" blocks.  */
#line 12 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.y"


#include "parsing/asm/common.h"
#include "parsing/asm/helper.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>


#line 68 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.tab.h"

/* Token kinds.  */
#ifndef ASM_YYTOKENTYPE
# define ASM_YYTOKENTYPE
  enum asm_yytokentype
  {
    ASM_YYEMPTY = -2,
    ASM_YYEOF = 0,                 /* "end of file"  */
    ASM_YYerror = 256,             /* error  */
    ASM_YYUNDEF = 257,             /* "invalid token"  */
    INST = 258,                    /* INST  */
    ID = 259,                      /* ID  */
    NUM = 260,                     /* NUM  */
    SEP = 261,                     /* SEP  */
    EOL = 262                      /* EOL  */
  };
  typedef enum asm_yytokentype asm_yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined ASM_YYSTYPE && ! defined ASM_YYSTYPE_IS_DECLARED
union ASM_YYSTYPE
{
#line 32 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.y"

    uint32_t                arg;
    char                    name[24];
    struct __asm_operand_t  operands;

#line 98 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.tab.h"

};
typedef union ASM_YYSTYPE ASM_YYSTYPE;
# define ASM_YYSTYPE_IS_TRIVIAL 1
# define ASM_YYSTYPE_IS_DECLARED 1
#endif


extern ASM_YYSTYPE asm_yylval;


int asm_yyparse (void);

/* "%code provides" blocks.  */
#line 23 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.y"


int yylex();

extern int  yyparse();
void        yyerror(char const*);


#line 122 "/Users/svittal/Documents/qontra/include/parsing/asm/parser.tab.h"

#endif /* !YY_ASM_YY_USERS_SVITTAL_DOCUMENTS_QONTRA_INCLUDE_PARSING_ASM_PARSER_TAB_H_INCLUDED  */
