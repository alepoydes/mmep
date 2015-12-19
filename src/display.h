#ifndef DISPLAY_H
#define DISPLAY_H

#include "vector.h"

extern real* display_buffer;
extern int is_aborting;
extern int is_new_frame;

int initDisplay(int* argc, char** argv);
void deinitDisplay();
void lockDisplay();
void releaseDisplay();
void displayRedraw();

// external keyboard handler
void keyboard_function(unsigned char key);
void mouse_function(int button, int state, real p[3]);
void motion_function(real p[3]);
void (*print_screen)(int width, int height, void* buffer);

#endif