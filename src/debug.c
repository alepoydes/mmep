#include <stdio.h>
#include <assert.h>

#include "debug.h"

void watch_number(real next, real prev, int digits) {
	char bufn[100],bufp[100]; 
	int nn=sprintf(bufn,"%.*"RF"f",digits,next); assert(nn<sizeof(bufn));
	int np=sprintf(bufp,"%.*"RF"f",digits,prev); assert(np<sizeof(bufp));
	//int ap=0; 
	int d=0;
	int c; for(c=0; c<nn && nn==np && bufn[c]==bufp[c]; c++) {
		if(bufn[c]>='0' && bufn[c]<='9') d++;
		//if(bufn[c]=='.') ap=1;
		fputc(bufn[c],stderr);
		if(d==DIGITS) fprintf(stderr,COLOR_FAINT);
	};
	if(prev>next) fprintf(stderr,COLOR_BLUE);
	else fprintf(stderr,COLOR_RED);
	for(; c<nn; c++) {
		if(bufn[c]>='0' && bufn[c]<='9') d++;
		fputc(bufn[c],stderr);
		if(d==DIGITS) fprintf(stderr,COLOR_FAINT);
	};
  	fprintf(stderr,COLOR_RESET);
}