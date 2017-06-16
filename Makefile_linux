CC = gcc

CFLAGS = -Wall -Werror \
	 -O3 -funroll-loops -march=native -mtune=native \
         -std=c99
GLUT = -lglut -lGLU -lGL

all: p7_interpolation

basic.o: basic.c basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
interpolation.o: interpolation.c interpolation.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 

p7_interpolation: p7_interpolation.o basic.o interpolation.o
	$(CC) $(CFLAGS) $^ -o $@ $(GLUT) -lm

p7_interpolation.o: p7_interpolation.c interpolation.h 
	$(CC) $(CFLAGS) -c $< -o $@


