#!/bin/bash
cp $PWD/exports/hg38_df38/hg38_df38-sequences.json $PWD/tests/test_exports/test_ex/test_ex-sequences.json
cp $PWD/exports/hg38_df38/hg38_df38-model_lengths.json $PWD/tests/test_exports/test_ex/test_ex-model_lengths.json
head -n 40000 $PWD/exports/hg38_df38/hg38_df38-mask.tsv > $PWD/tests/test_exports/test_ex/test_ex-mask.tsv
head -n 195260 $PWD/exports/hg38_df38/hg38_df38-byacc-full_region.tsv > $PWD/tests/test_exports/test_ex/test_ex-byacc-full_region.tsv
head -n 717 $PWD/exports/hg38_df38/hg38_df38-byacc-bench_region.tsv > $PWD/tests/test_exports/test_ex/test_ex-byacc-bench_region.tsv

cargo run -- --data-dir $PWD/tests/test_data/ --exp-dir $PWD/tests/test_exports/ prep-beds --assembly test_ex --in-tsv $PWD/tests/test_exports/test_ex/test_ex-mask.tsv --data-type masks
cargo run -- --data-dir $PWD/tests/test_data/ --exp-dir $PWD/tests/test_exports/ prep-beds --assembly test_ex --in-tsv $PWD/tests/test_exports/test_ex/test_ex-byacc-full_region.tsv --data-type assembly_alignments
cargo run -- --data-dir $PWD/tests/test_data/ --exp-dir $PWD/tests/test_exports/ build-idx --assembly test_ex --data-type assembly_alignments

mkdir $PWD/tests/test_data/test_ex/model_lengths
mkdir $PWD/tests/test_data/test_ex/sequences
cp $PWD/tests/test_exports/test_ex/test_ex-model_lengths.json $PWD/tests/test_data/test_ex/model_lengths/test_ex-model_lengths.json
cp $PWD/tests/test_exports/test_ex/test_ex-sequences.json $PWD/tests/test_data/test_ex/sequences/test_ex-sequences.json