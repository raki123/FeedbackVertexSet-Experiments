test_file="fast_tests"
PROGRAMNAME="./build/FeedbackVertexSet"
VERIFIER="./build/Verifier"

# CODE_NAME="\e[1m"
CODE_NAME="\033[1m"
CODE_OK="\e[32m"
CODE_WRONG="\e[31m"
CODE_TLE="\033[0;37m"
CODE_END="\033[0m"


#Just in case some test would loop.
MAX_TIME=900

cd build 
>&2 cmake -DCMAKE_BUILD_TYPE=Release ..
>&2 make 
cd ..


if [ $# -eq 0 ]; then
  >&2 echo "No arguments supplied, assuming test file = $test_file"
  >&2 echo "Possible usage ./collectTestsTimes [test_file]\n"
else
  test_file=$1
  >&2 echo "Reading tests from file $test_file"
  shift 1
fi

timetmp=`mktemp`
PROGRAM="$PROGRAMNAME $@"


while read p; do
  >&2 echo ${CODE_NAME}"Running $p" ${CODE_END}
  echo "=============================="
  echo "Test $p"
  /usr/bin/time -q --output=$timetmp -f "%E" timeout ${MAX_TIME}s $PROGRAM <../$p > ../${p%.graph}.out 2> a.err
  status=$?
  if [ $status -ne 0 ]; then
    echo ${CODE_NAME}${CODE_TLE}"[TIMEOUT]\n" ${CODE_END}
  fi
  printf "Time spent = "
  cat $timetmp
  $VERIFIER ../$p ../${p%.graph}.out ../${p%.graph}.opt
  echo "\n"
done <$test_file
