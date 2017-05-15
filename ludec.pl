#!usr/bin/perl
system('module load tools/mpich2-1.5-gcc');
#system('rm ludec.out ludec.err');
system('mpicc -o ludec ludec.c');
system('qsub ludec.pbs');

