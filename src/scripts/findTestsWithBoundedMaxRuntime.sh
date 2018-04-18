# Runs all tests and saves list of the tests that finish under MAX_TIME s

#cd ..
>&2 mkdir build
cd build
>&2 cmake -DCMAKE_BUILD_TYPE=Release ..
>&2 make 

FILES=../../tests/pace-16open/*.graph
OUTFILE="../fast_tests"
ile=0
MAX_TIME=2

if [ $# -eq 0 ]; then
  >&2 echo "No arguments supplied, assuming MAX_TIME = 2s"
  >&2 echo "Possible usage ./runAllTests [MAX_TIME]"
else
  MAX_TIME=$1
  >&2 echo "MAX_TIME set to $MAX_TIME seconds"
fi



rm $OUTFILE
touch $OUTFILE

totaltests=0
successtests=0
timetmp="mktemp"

for f in $FILES
do
  foo=${f#../../tests}
  >&2 echo "Processing $f file... "
  printf "$foo file...  ->  "
  /usr/bin/time --output=$timetmp -f "%E" timeout ${MAX_TIME}s ./FeedbackVertexSet <$f >a.out 2>a.err
  status=$?
 if [ $status -eq 0 ]; then
 	#echo "DOBRY Teścik"
 	cat $timetmp
 	successtests=$(( $successtests + 1 ))
 	foo=${f#../../}
 	echo $foo >> $OUTFILE
 else 
 	printf "TLE\n"
 fi
  #echo "Status = $status"
  totaltests=$(( $totaltests + 1 ))
done

echo "-----------------------"
echo "Passed $successtests out of $totaltests"
