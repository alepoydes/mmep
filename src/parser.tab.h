/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

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

#ifndef YY_YY_SRC_PARSER_TAB_H_INCLUDED
# define YY_YY_SRC_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NEG = 258,
    REAL = 259,
    SZ = 260,
    INTEGER = 261,
    FILENAME = 262,
    ID = 263,
    EOL = 264,
    SECS = 265,
    SECEF = 266,
    SECMA = 267,
    SECBC = 268,
    SECTV = 269,
    SECUC = 270,
    SECN = 271,
    SECEC = 272,
    SECDMV = 273,
    SECIMAGE = 274,
    SECLOADIMAGE = 275,
    SECPOSITIONS = 276,
    SECDIPOLE = 277,
    SECIMAGEFROMGNUPLOT = 278,
    SECCURRENT = 279,
    SECTEMP = 280,
    SECCUT = 281,
    SECDOMAIN = 282,
    TIP = 283,
    PLANE = 284,
    SPHERE = 285,
    UNIVERSE = 286,
    VERTEX = 287,
    UNIFORM = 288,
    RANDOM = 289,
    RELAX = 290,
    AFM = 291,
    BCFREE = 292,
    BCPERIODIC = 293
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 53 "src/parser.y" /* yacc.c:1909  */

 uint sz;
 int i;
 real r;
 real vec[3];
 char* fn;
 primitive_t* pr;

#line 102 "src/parser.tab.h" /* yacc.c:1909  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_SRC_PARSER_TAB_H_INCLUDED  */
