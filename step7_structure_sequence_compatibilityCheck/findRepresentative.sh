#!/bin/bash

# Check if a file is provided as an argument
if [ $# -ne 1 ]; then
	echo "Usage: $0 <input_file>"
	exit 1
fi

input_file="$1"
cluster=$(echo $input_file | sed -e 's/_hits_evaluated_gc15.tsv//g')
table=$(cat "${input_file}_hits_evaluated_gc15.tsv")

# Sort the table by column 9 (base pairs) in descending order
sorted_table=$(echo "$table" | grep -vi "seqID" | sort -k9,9gr | head -n30)

# Initialize variables to store the representative row and maximum values
representative_row=""
max_values=""

# Function to convert percentage values into a comparable string
compare_values() {
	local value1="$1"
	local value2="$2"
	
	# Convert comma-separated percentages to a comparable string
	IFS=',' read -ra v1 <<< "$value1"
	IFS=',' read -ra v2 <<< "$value2"

	for i in "${!v1[@]}"; do
			if (( $(echo "${v1[i]} > ${v2[i]}" | bc -l) )); then
				return 1
			elif (( $(echo "${v1[i]} < ${v2[i]}" | bc -l) )); then
				return 2
			fi
	done
	return 0
}

# Iterate through the sorted table and process column 8 (percentages)
while read -r line; do
	# Extract column 8 (percentage values)
	percentage_column=$(echo "$line" | awk '{print $8}')
	
	# Use IFS to split the percentage_column into individual percentages
	IFS=',' read -ra percentages <<< "$percentage_column"
	
	# Initialize a flag to check if all values are > 99%
	all_above_99=true

	for value in "${percentages[@]}"; do
		# Compare the number with 99 using bc
		if (( $(echo "$value < 99" | bc -l) )); then
			all_above_99=false
			break
		fi
	done

	# If all values are greater than 99%, check if this row has maximum values
	if [ "$all_above_99" = true ]; then
		if [ -z "$max_values" ]; then
			representative_row="$line"
			max_values="$percentage_column"
		else
			# Compare current row's percentage values with the current max values
			compare_values "$percentage_column" "$max_values"
			result=$?

			if [ $result -eq 1 ]; then
				representative_row="$line"
				max_values="$percentage_column"
			fi
		fi
	fi
done <<< "$sorted_table"

# Output the final representative row
if [ -n "$representative_row" ]; then
	echo -e "$cluster\t$representative_row";
else
	representative_row=$(echo "$sorted_table" | head -n1); #Assigning the top candidate having the best possible base-pairs for the model!
	echo -e "$cluster\t$representative_row";
fi
