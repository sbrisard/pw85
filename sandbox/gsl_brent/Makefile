all: brent_gsl.c
	gcc -c `gsl-config --cflags` -I/C:/opt/pw85/include -Wall brent_gsl.c
	gcc -o brent_gsl -L/mingw64/bin -LC:/opt/pw85/bin brent_gsl.o -lgsl -lm -lpw85
