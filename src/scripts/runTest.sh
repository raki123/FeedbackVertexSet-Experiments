#!/bin/bash

#cd ..
mkdir build
cd build
cmake .. 
make 

testName=$1

echo $testName
./FeedbackVertexSet < ../../tests/$testName