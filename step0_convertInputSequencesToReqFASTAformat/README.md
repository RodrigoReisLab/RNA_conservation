## Preparing input sequences for use in the pipeline

If you have the sequences in a multiline format as shown below.

>header1
AGTCGTAGCTACGTACGTACGTGTACGTACGTA
TGACGTACGTAGCTGCATGCTA

>header2
TGCAGATCGTAGTCGATGCTAGTGCATGCATGT
ACGTAGTGCAG

then, use the fasta_oneLiner.sh script to convert sequence in one-line.

Expected Output, below will be written to file ending in "_oneLine.fasta"
>header1
AGTCGTAGCTACGTACGTACGTGTACGTACGTATGACGTACGTAGCTGCATGCTA

>header2
TGCAGATCGTAGTCGATGCTAGTGCATGCATGTACGTAGTGCAG

Command on how to use the script is as follows 

    bash fasta_oneLiner.sh multiline.fasta

**Note:** If you have more than one file, you can use the fasta_oneLiner_job.sh script to execute the script.