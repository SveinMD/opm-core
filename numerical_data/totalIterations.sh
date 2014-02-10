#!/bin/bash
ext=".txt"
base1="testCase1_iterations-r-0"
base2="testCase1_iterations-r-00"
sumString=0

for i in {0..19}
do
    if [ 10 -gt $i ]
    then
        file=$base2$i$ext
    else
        file=$base1$i$ext
    fi
    sum=$(./sumIterations.py $file)
    sumString=$sumString+$sum
done

echo $sumString > sumstring.py
chmod 755 sumstring.py
./sumstring.py
