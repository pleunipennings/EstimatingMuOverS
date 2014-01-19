#!/bin/bash
mu=2e-04
cost=0.1
output_every_Xgen=50
numgen_inN=50.01
start_output=0.01
for seed in 1
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Code_and_shellscript/HIVevolution_HIV1site >./Code_and_shellscript/Link0.1_1.txt
done
