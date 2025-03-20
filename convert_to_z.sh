#!/usr/bin/env bash

work_dir="./work"
path_names="$work_dir/path_names.txt"
conv_plines="$work_dir/conv_plines.txt"
seq_out="$work_dir/seq_out.txt"
raw_zlines="$work_dir/raw_zlines.txt"
raw_zlines_path="$work_dir/raw_zlines_path.txt"
raw_qlines="$work_dir/raw_qlines.txt"
out_file="./$(basename $1).zlines.gfa"

mkdir -p $work_dir

echo "Converting p-lines from $1"
grep -E "^P" $1 | cut -f3 | sed -e "s/+/0/g" -e "s/-/1/g" -e "s/,/ /g" -e "s/$/ 0/" > $conv_plines
grep -E "^P" $1 | cut -f2 > $path_names

echo "Running sequitur"
cat $conv_plines | sequitur -d -e "0" -p > $seq_out

echo "Getting z-lines"
grep -E "^0 ->" $seq_out | sed -e "s/^0 -> //" -e 's/\[0\]/\n/g' -e "s/0\]/+]/g" -e "s/1\]/-]/g" > $raw_zlines_path
paste $path_names $raw_zlines_path | sed -e '/^\s*$/d' | sed -e "s/^/Z\t/" > $raw_zlines

echo "Getting q-lines"
grep -E "^[1-9]+ -> " $seq_out | sed -e "s/^/Q\t/" -e "s/ -> /\t/" > $raw_qlines

echo "Stitching things together"
grep -vE "^P" $1 > $out_file
cat $raw_zlines >> $out_file
cat $raw_qlines >> $out_file

