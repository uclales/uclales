#!/usr/bin/csh

mkdir ../runs/$fname
mv {$fname}[h.]* ../runs/$fname
cd ../runs
msrcp -R -pe 1000 $fname mss:CCBL
cd ../bin

