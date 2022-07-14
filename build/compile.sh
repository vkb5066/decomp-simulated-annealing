#!/usr/bin/csh

#Compile structure gen
mv Main_SA.c .tmp
gcc -o rgen -std=c11 -Ofast -w -lm *.c extern/kdtree-master/kdtree.c
mv .tmp Main_SA.c

#Compile simulated annealing
mv Main_SGen.c .tmp
gcc -o rsim -std=c11 -Ofast -w -lm *.c extern/kdtree-master/kdtree.c
mv .tmp Main_SGen.c
