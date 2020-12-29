#!/bin/bash

stage2_formating () {
    #formatiing the symbolic variable "VARij" to "VAR(i, j)"
    for ((i = 0; i < 10; i++))
    do
        for ((j = 0; j < 10; j++))
        do
            sed -i "s/$i$j/($i, $j)/g" $1
    	    #echo "($i, $j)"
        done
    done
}

stage2_formating "predict.stage1.txt"
stage2_formating "accelerometer_correct.stage1.txt"
stage2_formating "magnetometer_correct.stage1.txt"
