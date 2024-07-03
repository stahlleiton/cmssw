#!/bin/bash

# Check if the crab_dir argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <crab_dir>"
    exit 1
fi

# Get the crab_dir argument
crab_dir="$1"

# Check if the provided directory exists
if [ ! -d "${crab_dir}" ]; then
    echo "Error: Directory ${crab_dir} not found."
    exit 1
fi

# Loop through each subdirectory in crab_dir
for job_dir in ${crab_dir}/*; do
    # Check if the current path is a directory
    if [ -d "${job_dir}" ]; then
        # Get just the directory name without the path
        job_name=${crab_dir}$(basename "${job_dir}")
        # Run crab status for this job directory
		echo "====================================================================="
        echo "Querying status for job ${job_name}"
        crab status -d "${job_name}" --verboseErrors
        echo "====================================================================="
    fi
done