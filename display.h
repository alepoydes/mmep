#ifndef DISPLAY_H
#define DISPLAY_H

#include "vector.h"

extern real* display_buffer;
extern int is_aborting;

int initDisplay(int* argc, char** argv);
void deinitDisplay();
void lockDisplay();
void releaseDisplay();

#endif