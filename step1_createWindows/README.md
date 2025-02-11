## This is the first step, where the sequences if longer (e.g. genomes), will need to cut to obtain shorter sequences. 

The desired length of the sequences is user-defined.

We create windows using the script createWindows.pl. This script takes sequences in FASTA oneLine format and outputs windows in a file.

Command:

    perl createWindows.pl -f input_oneLine.fasta -w 100 -p 10 -O windows.fa

-w: length of the window.\
-p: user-defined length of overlap between two-consecutive windows. 
If -p set to 0, the windows will be non-overlapping.\
-O: output file name 


In the above example, input sequence will be cut into windows of length 100 with an overlap of 10nt.
**Note:** If you have more than one input file, you can use the windowing_job.sh script to execute the program.

**Important:** Merge all the windows file into 1 file.
	cat file1_windows.fasta file2_windows.fasta > all_windows.fasta

### Check and process windows for bad sequence content
In this step, we need to check if the windows obtained are correct, non-redundant and unique. Possible reasons can be
1. Due to sequencing errors, there may be sequences with Ns which has to be removed. Here, we discard any windows containing 'N' (ignoring case) in the sequence. 
2. Sequences with polyN (N=A|T|G|C) will not form base-pairs and hence can be discarded.
3. Windows can be identical (e.g. genomes from similar organisms, gene duplications, etc). These windows are redundant and thus only 1 is retained as a representative.

All this is taken care with 2 scripts: 1. removeNs_polyN_windows.pl and 2. removeDuplicates.pl. The scripts can be executed using process_windows_job.sh script. The output of script 1 is taken as input to script 2, resulting in a final output file ending with _uniq.fasta.
