#!/bin/sh
sed -e "s/catch ne_eo;/catch; ne_eo = lasterror;/" $1 > $1.sed
mv $1.sed $1

