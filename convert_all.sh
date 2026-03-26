for file in data/*.dna; do
    basename=$(basename $file .dna)
    sff parse $file -o ignore/$basename.json
    done