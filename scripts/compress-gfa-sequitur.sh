#!/usr/bin/env bash

work_dir="./work"
path_names="$work_dir/path_names.txt"
conv_plines="$work_dir/conv_plines.txt"
seq_out="$work_dir/seq_out.txt"
cleaned_seq_out="$work_dir/cleaned_seq_out.txt"
raw_zlines="$work_dir/raw_zlines.txt"
raw_zlines_path="$work_dir/raw_zlines_path.txt"
raw_zlines_path_fill="$work_dir/raw_zlines_path_fill.txt"
raw_qlines="$work_dir/raw_qlines.txt"
out_file="./$(basename $1).zlines.gfa"

mkdir -p $work_dir

echo "Converting p-lines from $1"
grep -E "^P" $1 | cut -f3 | sed -e "s/+/0/g" -e "s/-/1/g" -e "s/,/ /g" -e "s/$/ 0/" > $conv_plines
grep -E "^P" $1 | cut -f2 > $path_names

echo "Running sequitur"
cat $conv_plines | sequitur -d -e "0" -p -k "$2" > $seq_out

cat $seq_out | sed -e "s/\([^\[0-9]\)\([0-9]\)/\1>Q\2/g" -e "s/^\([0-9]\)/Q\1/" -e "s/\[\([0-9][0-9]*\)0\]/>\1/g" -e "s/\[\([0-9][0-9]*\)1\]/<\1/g" -e "s/ -> /\t/" -e "s/ //g" > $cleaned_seq_out

echo "Getting z-lines"
grep -P '^Q0\t' $cleaned_seq_out | sed -e "s/^Q0\t//" -e 's/\[0\]/\n/g' > $raw_zlines_path
cat $raw_zlines_path | sed -e '/^\s*$/d' -e "s/^/0\t0\t*\t*\t/" > $raw_zlines_path_fill
paste $path_names $raw_zlines_path_fill | sed -e '/^\s*$/d' | sed -e "s/^/Z\t/" > $raw_zlines

echo "Getting q-lines"
grep -P '^Q[1-9][0-9]*\t' $cleaned_seq_out | sed -e "s/^/Q\t/" > $raw_qlines

echo "Stitching things together"
grep -vE "^P" $1 > $out_file
cat $raw_zlines >> $out_file
cat $raw_qlines >> $out_file

