#!/bin/bash
# Script to complile and execute a C program

gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans
python3 setup.py build_ext --inplace