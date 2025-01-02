#!/bin/zsh
#
result_dir="../seq_align"
trimmed_dir="../queries/splitmixed"

for file in $trimmed_dir/*.fq; do
    echo "$file"
    sh sequence_alignment.sh $file &
done

wait
echo "All Done"
