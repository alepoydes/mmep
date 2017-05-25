/* recognize tokens for the calculator */
%option noyywrap nodefault yylineno
%{
# include "parser.h"
# include "parser.tab.h"
%}

%option nounput noinput noyywrap nodefault yylineno
%x IFILE

/* float exponent */
EXP ([Ee][-+]?[0-9]+)
%%

[\"] { BEGIN IFILE; }
<IFILE>[^\"]* { 
	yylval.fn=strdup(yytext);
	return FILENAME; 
}
<IFILE>\" { 
	BEGIN INITIAL;
}

"[size]" { return SECS; }
"[external field]" { return SECEF; }
"[magnetic anisotropy]" { return SECMA; }
"[boundary conditions]" { return SECBC; }
"[translation vectors]" { return SECTV; }
"[unit cell]" { return SECUC; }
"[neigborood]" { return SECN; }
"[exchange constant]" { return SECEC; }
"[dzyaloshinskii moriya vector]" { return SECDMV; }
"[image]" { return SECIMAGE; }
"[load image]" { return SECLOADIMAGE; }
"[load from gnuplot]" { return SECIMAGEFROMGNUPLOT; }
"[positions]" { return SECPOSITIONS; }
"[dipole]" { return SECDIPOLE; }
"[temperature]" { return SECTEMP; }
"[spin polarized current]" { return SECCURRENT; }
"[polarized current]" { return SECCURRENT; }
"[current]" { return SECCURRENT; }
"[cut]" { return SECCUT; }
"[domain]" { return SECDOMAIN; }

"relax" { return RELAX; }
"tip" { return TIP; }
"sphere" { return SPHERE; }
"plane" { return PLANE; }
"universe" { return UNIVERSE; }
"vertex" { return VERTEX; }
"uniform" { return UNIFORM; }
"random" { return RANDOM; }
"free" { return BCFREE; }
"periodic" { return BCPERIODIC; }
"afm" { return AFM; }

[A-Za-z][a-zA-Z0-9]* {
	yylval.fn=strdup(yytext);
	return ID;
}


"[" | 
"]" | 
"{" | 
"}" | 
"+" | 
"-" |
"!" | 
"*" |
"/" |
"(" |
")" |
"~" |
"," { return yytext[0]; }

-?[0-9]+"."[0-9]*{EXP}? |
"."[0-9]+{EXP}? { yylval.r = atof(yytext); return REAL; }
-?[0-9]+{EXP} { yylval.r = atof(yytext); return REAL; }

[0-9]+ { yylval.sz = atoi(yytext); return SZ; }

-[0-9]+ { yylval.i = atoi(yytext); return INTEGER; }

^[ \t\r]*"#".*"\n"
"#".*
^[ \t\r]*"\n"
[ \t\r] { /* ignore whitespace */ }

\n { return EOL; }

. { yyerror("Lexer error at line %d: Unexpected character '%c'\n", yylineno, *yytext); }

%%