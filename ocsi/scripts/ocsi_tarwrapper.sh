#!/bin/bash

#mark start
sleep 5
touch tempfile
sleep 5

#run thing
python /scripts/ocsi2.py "${@:1}"

#wrap up output and kick out tempfile
find . -mindepth 1 -newer tempfile -exec tar -rf FullOutput.tar {} \;
rm tempfile
