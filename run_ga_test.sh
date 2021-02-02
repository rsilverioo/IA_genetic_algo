
#/usr/bin/sh

prog=ga-test.exe
rm tests/*

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 # 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
	cp configs/GAconfig_$i.txt GAconfig
	echo Executing run $i of $prog with configuration file GAconfig_$i

	for j in 1 2 3 4 5 6 7 8 9 10
	do
			yes "" | ./$prog >> tests/results$i.txt
			echo "************************************" >> tests/results$i.txt

			mv raw_output.txt tests/raw_outputs/raw_output_{$i}_{$j}.txt
	done
done

# $SHELL