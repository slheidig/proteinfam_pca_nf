#!/bin/bash

# Script to collect all results for specific OGs into one folder
# Usage: ./collect_og_results.sh -i <input_dir> -o <output_dir> -t <target_OG>
# Example: ./collect_og_results.sh -i allOG -o OG334 -t OG0000334

# Default values
INPUT_DIR=""
OUTPUT_DIR=""
TARGET_OG=""

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) 
            INPUT_DIR="$2"
            shift
            ;;
        -o|--output) 
            OUTPUT_DIR="$2"
            shift
            ;;
        -t|--target) 
            TARGET_OG="$2"
            shift
            ;;
        *) 
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
    shift
done

# Validate arguments
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$TARGET_OG" ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 -i <input_dir> -o <output_dir> -t <target_OG>"
    echo "Example: $0 -i allOG -o OG334 -t OG0000334"
    exit 1
fi

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Collecting files for ${TARGET_OG}..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""

# Find all files matching the OG pattern in the input directory
# This includes files with OG in filename or in directory path
find "$INPUT_DIR" -type f \( -name "*${TARGET_OG}*" -o -path "*${TARGET_OG}/*" \) | while read -r file; do
    # Get just the filename
    filename=$(basename "$file")
    dest="$OUTPUT_DIR/$filename"
    
    # Check if file already exists in output directory
    if [[ -e "$dest" ]]; then
        # Add counter suffix to avoid overwriting
        counter=1
        # Handle files with and without extensions
        if [[ "$filename" == *.* ]]; then
            name="${filename%.*}"
            ext="${filename##*.}"
            while [[ -e "$OUTPUT_DIR/${name}_${counter}.${ext}" ]]; do
                ((counter++))
            done
            dest="$OUTPUT_DIR/${name}_${counter}.${ext}"
        else
            while [[ -e "$OUTPUT_DIR/${filename}_${counter}" ]]; do
                ((counter++))
            done
            dest="$OUTPUT_DIR/${filename}_${counter}"
        fi
        echo "  Warning: File '$filename' already exists, saving as '$(basename "$dest")' to avoid overwriting"
    fi
    
    # Copy file to output directory
    cp "$file" "$dest"
    echo "  Copied: $(basename "$dest")"
done

echo ""
echo "All files for ${TARGET_OG} collected in: $OUTPUT_DIR"
