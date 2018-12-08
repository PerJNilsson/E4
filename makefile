.PHONY : all

all : compile

CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

compile :
	$(CC) main.c -o run $(LIBS) $(CFLAGS)

plot :
	python plot_velocities.py; python plot_positions.py; python plot_hist_x.py; python plot_hist_v.py;

clean :
	rm -rf run; rm -rf *.dat;
