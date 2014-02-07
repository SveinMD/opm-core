#!/bin/bash
file="testCase"
ext=".out"
outfile=$file$ext
sum=0
for i in {1..10}
do
 #outfile=$file$i$ext
 ./testCase1 -s t &> $outfile
 resultFull=$(tail -1 $outfile | head -1)
 result=$(echo $resultFull java | cut -d ' ' -f4)
 sum=$sum+$result
done
echo $sum
