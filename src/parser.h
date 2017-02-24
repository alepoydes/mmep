#ifndef PARSE_H
#define PARSE_H

#include "vector.h"
#include "cut.h"

#include <stdio.h>

/* interface to the lexer */
extern int yylineno; /* from lexer */
extern int yydebug;
void yyerror(const char *s, ...);

int yylex(void);

void parse_lattice(FILE* file);

#endif