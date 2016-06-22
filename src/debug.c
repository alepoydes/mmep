#include <stdio.h>
#include <assert.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>

#include "debug.h"

volatile int stop_signal=0;

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
};

#define MAX_SIG 3

void signal_handler(int sig) {
	if(sig==SIGINT) {
    	fprintf(stderr, COLOR_YELLOW"received SIGINT %d/%d\n"COLOR_RESET, ++stop_signal, MAX_SIG);
    	if(stop_signal>MAX_SIG) exit(1);
	};
};

void init_signal() {
	if(signal(SIGINT, signal_handler)==SIG_ERR)
		fprintf(stderr, COLOR_RED COLOR_BOLD"can't catch SIGINT\n"COLOR_RESET);
};