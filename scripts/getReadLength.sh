in_file=$1

first_line=`zcat $1 | head -n 2 | tail -n 1`
read_len=`expr length "${first_line}"`

echo $read_len