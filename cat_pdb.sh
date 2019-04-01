#!/bin/sh

for f in data/*.pdb;do 
   cat $f | sort -n -k 4 >> all.pdb
   echo END >> all.pdb
   echo $f
done
mv -v all.pdb output.pdb




