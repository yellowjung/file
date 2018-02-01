#!/bin/bash

autoreconf -f -i
CC=clang CFLAGS='-O3 -mtune=haswell -flto' LDFLAGS='-flto -lm' ./configure --disable-zlib --prefix=/tmp/file && make -j$(grep -c ^processor /proc/cpuinfo)
