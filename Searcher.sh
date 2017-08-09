#!/bin/bash
echo "Enter input file name...";
read fname;
echo "Enter the search pattern";
read pattern
for x in $(cat $fname);
do
	if [ -f $x ]; then
		for y in $(cat $pattern);
		do
			echo "$y" >> ${x}_counts.txt
			`grep -c -w "$y" $x >> ${x}_counts.txt`
		done
	fi
done
