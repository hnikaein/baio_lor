#!/usr/bin/env bash
output_file_name=outputs/te26
read_file_name=grch38.reads3.fastq

rm $output_file_name
rm c_long_read_aligner
cmake --build cmake-build-debug/ --target c_long_read_aligner -- -j 8 && mv cmake-build-debug/c_long_read_aligner .
chmod +x c_long_read_aligner
./c_long_read_aligner -r seqs/grch38.fasta -q seqs/$read_file_name -l 5 -o $output_file_name -s &