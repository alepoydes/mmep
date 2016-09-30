#ifndef CMD_H
#define CMD_H

#include <stdio.h>

extern real epsilon;
extern long int max_iter;
extern real mode_param;
extern real param2;
extern int mode;
extern int integrator;
extern int debug_plot;
extern int debug_every;
extern int save_octave;
extern real dipole_negligible;
extern const char* outdir;

void print_settings();

int init_program(int argc, char** argv,
const char* info, const char* options_desc,
const char* optstr, char (*handle)(char opt, const char* arg)
);

void oct_save_lattice(FILE* file);

#endif