>&2 mkdir build
cd build 
>&2 cmake -DCMAKE_BUILD_TYPE=Release ..
>&2 make Verifier
cd ..


FILES=../tests/pace-16open/*.opt
PROGRAM="./build/FeedbackVertexSet"
VERIFIER="./build/Verifier"

# CODE_NAME="\e[1m"
CODE_NAME="\033[1m"
CODE_OK="\033[32m"
CODE_WRONG="\033[31m"
CODE_TLE="\033[0;37m"
CODE_END="\033[0m"


failedtests=0
for f in $FILES
do
  >&2 echo ${CODE_NAME}"Running $f" ${CODE_END}
  foo=${f#../}
  echo "=============================="
  echo "Test $foo"

  echo "$VERIFIER ../${foo%.opt}.graph ../${foo%.opt}.out ../$foo"
  $VERIFIER ../${foo%.opt}.graph ../${foo%.opt}.out ../$foo
  status=$?

  if [ $status -ne 0 ]; then
    failedtests=$(( $failedtests + 1 ))
  fi

  # $VERIFIER ../${p%.out}.graph ../${p%.graph}.out ../$p
  # echo "\n"
done

echo "\n--------------------------------------------\n"

if [ $failedtests -ne 0 ]; then
  echo ${CODE_NAME}${CODE_WRONG}"Error in $failedtests output files"${CODE_END}
else
  echo ${CODE_NAME}${CODE_OK}"All outputs are correct" ${CODE_END}
fi
