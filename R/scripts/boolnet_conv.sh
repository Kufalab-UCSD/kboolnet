# Small script to make boolnet files human readable by replacing node IDs with rxncon names
# Usage: boolnet_conv.sh <_symbols.csv> <.boolnet>

# Generate regexes from the symbols file
regexes="$(sed 's|\(.*\), \(.*\)|s/\1/\2/g|' $1)"

# Split regexes into array 
SAVESAVEIFS=$IFS   # Save current IFS
IFS=$'\n'      # Change IFS to new line
regexes=($regexes) # split to array $names
IFS=$SAVEIFS   # Restore IFS

# Read boolnet file into variable
boolnet="$(cat $2)"

# Apply regexes
for exp in "${regexes[@]}"
do
    boolnet="$(echo $boolnet | sed $exp)"
done

# Output new boolnet
echo $boolnet
