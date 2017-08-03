#!/bin/sh
#
# replace the "catch EXPRESSION;" with plain contstruct

find ../.. -name "*.m" -type f -exec ./replace_catch_ne_eo.sh {} \;
find ../.. -name "*.?ff" -type f -exec ./replace_catch_ne_eo.sh {} \;
