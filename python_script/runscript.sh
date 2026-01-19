#!/bin/bash

test_run()
{
  paramfolder="$1"
  cd "$paramfolder"
  echo python ../main.py "$paramfolder.json"
  cd ..
}
run()
{
  paramfolder="$1"
  cd "$paramfolder"
  python ../main.py "$paramfolder.json"
  cd ..
}
if [ -d "$1" ]; then
	test_run "$1"
	echo "Enter to continue, Ctrl+c to end"
	read option
	if [[ -z $option ]]; then
	 	run "$1"
	else
		echo "run program denied!"
	fi 
else
    echo "give argv as param folder name"
fi