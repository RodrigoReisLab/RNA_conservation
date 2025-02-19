# This step is to further validate the compatibility of sequences with the structure.

This step is performed using nhmmer (It is already present along with R-scape). When we perform RNA structure prediction using LocARNA, it may occur that the sequence alignments are forced to occupy positions which contain structure and thus, the alignments may not be accurate, giving rise to false covariation in base-pairs. Thus, in this step we correct the alignments, remove sequences which do not fit the alignment and obtain a better alignment representing the structure. 

The alignment obtained from Step 6 is taken, it is first de-aligned, removed with structure information carried out using esl-reformat

	~/Documents/rscape_v2.0.4.a/lib/hmmer/easel/miniapps/esl-reformat -o motif.fa fasta motif.sto

Here, the input file motif.sto is obtained from previous step (can either be motif_cleaned.sto or motif_cleaned.cacofold.sto). The output will be written as fasta unaligned, ungapped sequences.  

A representative sequence is then chosen from the sequences. Here, the output generated from RNA-SCoRE (motif_evaluated.tsv) in Step 6.iii will guide you to choose the best candidate. Ideally, you want to choose the candidate sequence which has passed the evaluation, is the longest, has maximum base-pairing probability and maximum motif present. The script 'findRepresentative.sh' can be used as a guide to find the best candidate. However, you can also select one manually, since it is not the most efficient way to find representative. Also, since there is no right or wrong selection, the script-chosen candidate may not be the one you think is the best candidate.  

This output should written to a tab-separated file 'listOfRepresentative', where column-1 indicates the name of the motif and column-2 indicates header of the representativesequence (from FASTA file).

<pre>
motif_cleaned			NC_040945.1/52289-52360_1
NC_036423.1_1_175_425r		NC_036980.1/437-565_-1
NC_036093.1_258_45150_45400r	NC_033344.1/176-235_-1
NC_029965.1_790_138250_138500r	NC_043954.1/141445-141658_1
NC_028519.1_387_67725_67975r	NC_009267.1/67748-67916_-1
NC_035414.1_158_27650_27900r	NC_022408.1/32753-32915_-1
NC_036113.1_442_77350_77600r	NC_030621.1/15347-15520_1
NC_043805.1_34_5950_6200r	NC_029961.1/5905-6060_-1
</pre>

Once the input files are ready, we perform several tasks:
i. Get the sequence of the representative in a separate FASTA file and convert to an alignment stockholm file. This will act as a model to align rest of the sequences.  
Representative .sto file should like this. Sequence and SS_cons should be a single line. 

<pre>
# STOCKHOLM 1.0
#=GF ID NC_042825.1_109_19075_19325r
#=GF AU Infernal 1.1.4

NC_040945.1/52289-52360_1 CACCUACUUAACUCAGCGGUUAGAGUAUCGCUUUCAUACGGCGGGAGUCAUUGGUUCAAAUCCAAUAGUAGG
#=GC SS_cons              ::((((((,,<<<<________>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))
//
</pre>

ii. Build model from the representative using hmmbuild and then use the model to search against rest of the sequences using hmmsearch.  

iii. Find covariation just based on sequences and use --cacofold to obtain structure for the newly obtained alignment.  

The script 'nhmmer_local_job.sh' can then be used as a guide to perform all the above steps. The output will be a .sto alignment file with R-scape cacofold predicted structure. Lets call this as _nhmmer_Rscape_cacfold_structure_.

As a sanity check, add the locarna predicted structure to the alignment. Since the sequences were already from a structure based alignment, all the positions in the hmmer-derived alignment with gaps are added with gaps in the LocARNA predicted structure. This is performed using the script 'addStructToHMMERaln.pl'. Make sure the representative .sto file is not a multi-line (wrapped) sequence but in a single line. Else the program will not assign the structure correctly.  

After adding the structure perform R-scape with two-tailed test (_-s_ option). Lets call this as _locarnaAddedStructureToHmmAln_. Compare the structures _nhmmer_Rscape_cacfold_structure_ and _locarnaAddedStructureToHmmAln_. Use the same logic as previously described in Step 6.iv and you will obtain reliable alignment compatible with structure. Additionally, if the covariation is lost in this process, the structures may become less reliable and hence should be discarded.

###This way, with all the steps scrutinizing the sequences, the structure and the compatibility between sequence and structure, the possibility of obtaining reliable models increases and helps in predicting high-confidence de novo motifs in the sequences of interest.


### Steps after selecting motifs:
Ideally, once the motifs are obtained, these can be used as an improved seed to expand your search for homologs in various other sequences/ genomes of interest (Scaling up). If the covariation does not improve from the previous step or you find that this model is performing worse than what you obtained originally, you can use the same steps on different set of sequences as database to find new homologs to improve seed. Thus, iteratively performing steps 5-7 (see Figure 1: flowchart) until the user is satisfied with the alignment, the model and the homologs. One approach to test if the models are high-confidence, is if no new hit is obtained even after several iterations.

The best alignment can then be considered seed and all the homologs identified by the seed is called as family. These can be submitted to [Rfam](https://rfam.org/).  
