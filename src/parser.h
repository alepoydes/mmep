#ifndef PARSE_H
#define PARSE_H

#include <stdio.h>
#include "vector.h"

/* interface to the lexer */
extern int yylineno; /* from lexer */
extern int yydebug;
void yyerror(const char *s, ...);

int yylex(void);

void parse_lattice(FILE* file);

#endif