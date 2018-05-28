#!/bin/bash
file="dynmat.mold";
flag1="FREQ"
flag2="FR-COORD"
grep -A 100 $flag1 $file | grep -B 100 $flag2 | grep -Eo "[0-9]{0,9}[\.]{1}[0-9]{0,9}"
