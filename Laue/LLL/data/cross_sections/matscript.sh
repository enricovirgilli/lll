#!/bin/bash
# bash shell script

for mat in CdTe GaAs
do
pyxcom --range [0000:1000] --mat $mat --columns [1,7,8]|grep -v "#"|sed "s/inf/0/g" > $mat.dat
pyxcom --range [1000:2000] --mat $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [2000:3000] --mat $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [3000:4000] --mat $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [4000:5000] --mat $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
done
