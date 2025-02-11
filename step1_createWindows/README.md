## This is the first step, where the sequences if longer (e.g. genomes), will need to cut to obtain shorter sequences. 

The desired length of the sequences is user-defined.

We create windows using the script createWindows.pl. This script takes sequences in FASTA oneLine format and outputs windows in a file.

Command:

    perl createWindows.pl -f input_oneLine.fasta -w 100 -p 10 -O windows.fa

-w: length of the window.\
-p: user-defined length of overlap between two-consecutive windows. 
If -p set to 0, the windows will be non-overlapping.\
-O: output file name 

**Note:** If you have more than one input file, you can use the windowing_job.sh script to execute the program.

### Check and process windows for bad sequence content

