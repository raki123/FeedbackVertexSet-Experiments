test_file="fast_tests"
PROGRAM="./build/Stats"


#Just in case some test would loop.
MAX_TIME=900

cd build 
>&2 cmake -DCMAKE_BUILD_TYPE=Release ..
>&2 make 
cd ..


if [ $# -eq 0 ]; then
  >&2 echo "No arguments supplied, assuming test file = $test_file"
  >&2 echo "Possible usage ./statsCollect  [test_file]\n"
else
  test_file=$1
  >&2 echo "Reading tests from file $test_file"
  shift 1
fi

timetmp=`mktemp`


while read p; do
  /usr/bin/time -q --output=$timetmp -f "%E" timeout ${MAX_TIME}s $PROGRAM <../$p > ../${p%.graph}.stats 2> a.err
  stats=`head -n 1 ../${p%.graph}.stats | cut -f1`
  echo "$p,$stats"
done <$test_file
