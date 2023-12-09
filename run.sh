#!/bin/bash

command_to_run="python ./benchmark/profile.py space_complexity_info_updated.txt benchmarking_results_updated.csv > space_complexity_info_updated.txt"

# Run the command 25 times
for ((i=1; i<=25; i++)); do
    echo "Running command iteration $i"
    $command_to_run
done