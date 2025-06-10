#!/bin/bash

# Loop through each file in the current directory
for file in *; do
    # Skip directories and the script itself
    if [[ -d "$file" || "$file" == "$0" ]]; then
        continue
    fi
    
    
    # Check if the file is a Java file
    if [[ "$file" == *.zip ]]; then
        # Extract the part after the last underscore
        
        unzip "$file" -d "data/."
    fi

echo "Files have been placed inside of data dir."
done