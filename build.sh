#!/bin/bash

gcc -O3 -c pileup_stat_lib.c   -W -Wall
gcc -O3 -o count count.c pileup_stat_lib.o  -lm -W -Wall
