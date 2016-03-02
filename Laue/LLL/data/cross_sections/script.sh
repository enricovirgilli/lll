#!/bin/bash
# bash shell script

for mat in `seq 99`
do
pyxcom --range [0000:1000] -Z $mat --columns [1,7,8]|grep -v "#"|sed "s/inf/0/g" > $mat.dat
pyxcom --range [1000:2000] -Z $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [2000:3000] -Z $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [3000:4000] -Z $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
pyxcom --range [4000:5000] -Z $mat --columns [1,7,8]|grep -v "#" >> $mat.dat
done