#!/bin/bash

# Usage : ./iter.sh /path/to/iterate > /path/to/save/log
# Sort  : cat log | grep -v empty | sort -k3 -g

FILE=./src/.libs/file
MAGIC=./magic/magic.mgc
find "$@" -xdev -type f | xargs -n8192 -P$(grep -c ^processor /proc/cpuinfo) -d '\n' $FILE -m $MAGIC -f -
