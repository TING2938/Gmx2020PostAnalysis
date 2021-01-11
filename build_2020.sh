#!/usr/bin/env bash

## usage
# source path/to/GMXRC
# ./build_2020.sh path/to/your/code

basename=${1%%.c*}
g++ $* -o $basename -I $GMXPREFIX/include -I ./include -L $GMXLDLIB -lgromacs -O3 -std=c++17

