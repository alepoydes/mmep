#ifndef DEBUG_H
#define DEBUG_H

#include "stdio.h"
#include "vector.h"

#define COLOR_RESET "\033[m"
#define COLOR_BOLD "\033[1m"
#define COLOR_FAINT "\033[2m"
#define COLOR_ITALIC "\033[3m"
#define COLOR_RED "\033[31m"
#define COLOR_GREEN "\033[32m"
#define COLOR_YELLOW "\033[33m"
#define COLOR_BLUE "\033[34m"
#define CLEAR_SCREEN "\033[2J"
#define RESET_CURSOR "\033[1;1H"

extern volatile int stop_signal;

#define TRUE 1
#define FALSE 0

FILE* open_file(const char* dir, const char* filename, char to_write);

void print_vector(int N, real* data);
void watch_number(realp next, realp prev, int digits);
void fprint_timediff(FILE* file, double timediff);

void init_signal();

#endif