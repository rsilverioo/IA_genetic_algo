
#/usr/bin/sh

prog=ga-test-report.exe

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	yes "" | ./$prog $i
	echo Executing run $i of $prog with configuration file GAconfig
done

# $SHELL