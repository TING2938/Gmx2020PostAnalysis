#!/usr/bin/env bash

source path/to/GMXRC

basename=${1%%.c*}
g++ $* -o $basename -I $GMXPREFIX/include -L $GMXLDLIB -lgromacs -O3 -std=c++17

