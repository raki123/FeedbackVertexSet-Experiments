test_file="fast_tests"
PROGRAM="./build/FeedbackVertexSet Appx"
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


while read p; do
  /usr/bin/time -q --output=$timetmp -f "%E" timeout ${MAX_TIME}s $PROGRAM <../$p > ../${p%.graph}.appx 2> a.err
  status=$?
  >&2 $VERIFIER ../$p ../${p%.graph}.appx ../${p%.graph}.opt
  apx=`head -n 1 ../${p%.graph}.appx | cut -f1`
  if [ $status -ne 0 ]; then
    t=`RTE`
  else
    t=`cat $timetmp | cut -f1`
  fi
  echo "$p,$apx,$t"
done <$test_file
