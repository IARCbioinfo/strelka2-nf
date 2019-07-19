#!/bin/bash

for vcf in "$@"
do

out=${vcf/vcf.gz/vcf}

first_format_num=$(zgrep -n -m 1 '##FORMAT' "$vcf" | cut -d : -f 1)
zcat "$vcf" | sed "$first_format_num"'i##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' > "$out"
sed -ri 's|(DP:)|GT:\1|g' "$out"
sed -ri 's|(:TU\t)|\10/0:|g' "$out"
sed -ri 's|(:TU\t[^\t]*\t)|\10/1:|g' "$out"
sed -ri 's|(:BCN50\t)|\10/0:|g' "$out"
sed -ri 's|(:BCN50\t[^\t]*\t)|\10/1:|g' "$out"
gzip -f $out

done
