all: min min2
#OPT=-std=c99 -flto -Ofast -march=native
OPT=-std=c99 -g -Wall

min: parse.c vector.c skyrmion.c plot.c optim.c min.c 
	gcc $^ $(OPT) -lm -o $@	

min2: parse.c vector.c skyrmion.c plot.c optim.c min2.c 
	gcc $^ $(OPT) -lm -o $@	
