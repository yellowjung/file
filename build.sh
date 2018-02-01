#!/bin/bash

# sed -i -e 's@# deb-src@deb-src@g'
# apt update
# apt build-dep file

autoreconf -f -i
CC=clang CFLAGS='-O3 -mtune=haswell -flto' LDFLAGS='-flto -lm' ./configure --disable-zlib --prefix=/tmp/file && make -j$(grep -c ^processor /proc/cpuinfo)
