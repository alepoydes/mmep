#include <stdio.h>
#include <assert.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "debug.h"

volatile int stop_signal=0;

FILE* open_file(const char* dir, const char* filename, char to_write) {
	char buf[1024]="";
	if(dir) strncat(buf, dir, sizeof(buf)-strlen(buf)-1);
	if(filename) strncat(buf, filename, sizeof(buf)-strlen(buf)-1);
    fprintf(stderr, COLOR_YELLOW "%s" COLOR_RESET " '%s'\n",to_write?"Saving":"Reading",buf);
	FILE* file=fopen(buf,to_write?"wb":"rb");
	if(!file) fprintf(stderr, COLOR_RED "Can not open '%s' for %s\n" COLOR_RESET,buf,to_write?"writing":"reading");
	return file;
};

void fprint_timediff(FILE* file, double timediff) {
	int f=0;
	int sec=timediff;
	int min=sec/60; sec%=60;
	int hour=min/60; min%=60;
	int day=hour/24; hour%=24;
	if(day!=0) { fprintf(file, "%dd", day); f=3; };
	if(f>0 || hour!=0) { fprintf(file, "%dh", hour); f=2; };
	if(f>0 || min!=0) { fprintf(file, "%dm", min); f=1; };
	if(f<2) fprintf(file, "%ds", sec); 
};

void watch_number(realp next, realp prev, int digits) {
	char bufn[100], bufp[100]; 
	uint nn=sprintf(bufn,"%.*" RPF "f",digits,next); assert(nn<sizeof(bufn));
	uint np=sprintf(bufp,"%.*" RPF "f",digits,prev); assert(np<sizeof(bufp));
	//int ap=0; 
	uint d=0;
	uint c; for(c=0; c<nn && nn==np && bufn[c]==bufp[c]; c++) {
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

#define MAX_SIG 2

void signal_handler(int sig) {
	if(sig==SIGINT) {
    	fprintf(stderr, COLOR_YELLOW "received SIGINT %d/%d\n" COLOR_RESET, ++stop_signal, MAX_SIG);
    	if(stop_signal>=MAX_SIG) exit(1);
	};
};

void init_signal() {
	if(signal(SIGINT, signal_handler)==SIG_ERR)
		fprintf(stderr, COLOR_RED COLOR_BOLD "can't catch SIGINT\n" COLOR_RESET);
};

void print_vector(int N, real* data) {
	printf("[");
	for(int n=0; n<N; n++) 
		printf("%s%" RF "g", n>0?COLOR_FAINT "," COLOR_RESET:"", data[n]);
	printf("]\n");
};