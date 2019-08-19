#!/bin/bash
# Remember to chmod u+x
# use nohup
while ps -p 4140 >/dev/null 2>&1
do
   echo sleeping
   sleep 30m
done 
./Release/SqWWW /Users/alggroup/Documents/SqWWW/Inputs/lvt_005.txt
./Release/SqWWW /Users/alggroup/Documents/SqWWW/Inputs/mot_005.txt
wait
./Release/SqWWW /Users/alggroup/Documents/SqWWW/Inputs/qzd_005.txt
./Release/SqWWW /Users/alggroup/Documents/SqWWW/Inputs/usf_005.txt
