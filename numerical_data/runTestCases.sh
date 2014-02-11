#!/bin/bash

float_scale=2

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0 && -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat 
}

file="testCase"
ext=".out"
outfile=$file$ext
e="0"

echo "Solver [r/t/i/b/n + return]:"
read solver

echo "Print [-p/<blank> + return]:"
read do_print

echo "**Running solvers**"
echo "..."

for i in {1..10}
do
 #outfile=$file$i$ext
 #../bin/testCase1 -p -s b &> $outfile
 ../bin/testCase1 $do_print -s $solver &> $outfile
 resultFull=$(tail -1 $outfile | head -1)
 result=$(echo $resultFull java | cut -d ' ' -f4)
 e=$e" + "$result
done
echo "Test run in $(float_eval "$e") seconds"
