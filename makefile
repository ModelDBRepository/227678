all: nmodel fixed clean

atp.o : atp.h atp.c
	gcc  -c atp.c

euler_noise.o: atp.h euler_noise.c atp.c
	gcc -c euler_noise.c -lm -Wmultichar

nmodel:  euler_noise.o
	gcc  -W -Wall -Wmultichar -o nmodel euler_noise.o -lm

fixed_finder.o: atp.h fixed_finder.c
	gcc -c fixed_finder.c -lm

fixed:  fixed_finder.o fixed_finder.c atp.h
	gcc -W -Wall -o fixed fixed_finder.o -lm

clean:
	$(RM) *.o
