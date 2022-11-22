
#Anellovirus_2022.0.fasta #$1

echo "global db" > Dicodb.py
echo -n "db = { " >> Dicodb.py
for headerInfo in $(grep ">" Anellovirus_2022.0.fasta | cut -d">" -f2):
 do echo -n "'"$(echo $headerInfo | cut -f1 -d",")"'"": [ '"$(echo $headerInfo | cut -f2 -d",")"','"$(echo $headerInfo | cut -f3 -d",")"','"$(echo $headerInfo | cut -f4 -d",")"','"$(echo $headerInfo | cut -f5 -d",")"','"$(echo $headerInfo | cut -f6 -d",") "','" >> Dicodb.py; done
 for seq in $(grep -A1 $(echo $headerInfo | cut -f1 -d",") Anellovirus_2022.0.fasta) ; do if [ ${seq::1} != ">" ] ; then echo -n $seq"']," >> Dicodb.py ; fi ;  done
echo "'none': [] }" >> Dicodb.py 



