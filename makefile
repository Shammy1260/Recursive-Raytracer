raytracer:	raytracer.o
		g++ -std=c++11	-Wall	-O	raytracer.o	-o	raytracer

raytracer.o:	raytracer.cpp	raytracer.h
				g++	-std=c++11 -c	-O	raytracer.cpp

clean:
			rm *.o *.ppm raytracer
