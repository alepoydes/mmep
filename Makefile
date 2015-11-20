all: min mind minq min2 min2d min2q mep mepd mepq 
OPT=-std=c99 -flto -Ofast -march=native -fopenmp
#OPT=-std=c99 -flto -Ofast -march=native
#OPT=-std=c99 -g -Wall -DDEBUG 

min: parse.c vector.c skyrmion.c plot.c optim.c min.c 
	gcc $^ $(OPT) -lm -o $@	

mind: parse.c vector.c skyrmion.c plot.c optim.c min.c 
	gcc $^ $(OPT) -DDOUBLE -lm -o $@	

minq: parse.c vector.c skyrmion.c plot.c optim.c min.c 
	gcc $^ $(OPT) -DQUAD -lm -o $@	

min2: parse.c vector.c skyrmion.c plot.c optim.c min2.c 
	gcc $^ $(OPT) -lm -o $@	

min2d: parse.c vector.c skyrmion.c plot.c optim.c min2.c 
	gcc $^ $(OPT) -DDOUBLE -lm -o $@	

min2q: parse.c vector.c skyrmion.c plot.c optim.c min2.c 
	gcc $^ $(OPT) -DQUAD -lm -o $@	

mep: parse.c vector.c skyrmion.c plot.c optim.c mep.c 
	gcc $^ $(OPT) -lm -o $@	

mepd: parse.c vector.c skyrmion.c plot.c optim.c mep.c 
	gcc $^ $(OPT) -DDOUBLE -lm -o $@	

mepq: parse.c vector.c skyrmion.c plot.c optim.c mep.c 
	gcc $^ $(OPT) -DQUAD -lm -o $@	