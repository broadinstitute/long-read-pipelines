if [[ $1 =~ \.gz$ ]]
then gzcat $1 | grep -v ^# | wc -l
else cat $1 | grep -v ^# | wc -l
fi