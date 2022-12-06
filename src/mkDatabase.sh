#!/bin/bash
#$1 path to database

echo "global db" > Dicodb.py
echo -n "db = { " >> Dicodb.py
for lane in $(cat $1):
  do if [ ${lane::1} != ">" ] ;
     then echo -n $lane"']," >> Dicodb.py;
     else headerInfo=$(echo $lane | cut -d">" -f2);echo -n "'"$(echo $headerInfo | cut -f1 -d",")"'"": [ '"$(echo $headerInfo | cut -f2 -d",")"','"$(echo $headerInfo | cut -f3 -d",")"','"$(echo $headerInfo | cut -f4 -d",")"','"$(echo $headerInfo | cut -f5 -d",")"','"$(echo $headerInfo | cut -f6 -d",") "','" >> Dicodb.py;
  fi;
done;
echo "'none': [] }" >> Dicodb.py