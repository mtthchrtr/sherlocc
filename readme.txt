#----------------------------------------------#
#         -----   SHERLOCC 1.0   -----         #
#----------------------------------------------#

SHERLOCC - A simple tool to detect conserved rare codon clusters in protein families and predict family consensus unstructured regions.

Designed by Matthieu Chartier and Rafael Najmanovich. Programmed by Matthieu Chartier.

Najmanovich Research Group - Faculty of medicine and health sciences, Biochemistry Departement, Universite de Sherbrooke, Sherbrooke, Quebec, Canada.
bcb.med.usherbrooke.ca

Sherlocc allows any user to search Pfam protein family alignments for regions with conserved low codon usage frequency. It also uses the equation of Uversky [1] with the window averaging method proposed by Prilusky [2] to predict unstructured regions in the protein family.



#------------------------------------------------#
#------How to analyse a Pfam protein family------#
#------------------------------------------------#

Sherlocc lets you use input Pfam alignment files that you can directly download from the Pfam website or alternatively you can start the analysis from Stage 3 using as input the precalculated Stage2 output files available on our website. Starting from stage 3 is much faster since all the web queries are skipped (which is the limiting step). Upon running sherlocc.pl from the terminal you will be prompted to start from Stage 3 or from Stage 1.

You can Download Stage2 output files from the BCB website (http://bcb.med.usherbrooke.ca).


-------------------------------------
@@----If you start from stage 1----@@
-------------------------------------

1. Download a Pfam family alignment from http://pfam.sanger.ac.uk/. The seed and the full alignments can be used. Since Sherlocc does web queries to retrieve information about the proteins of the alignment, full alignment that contain many proteins can take some time to analyse.

2. Place the uncompressed Pfam alignment file in pfam_files/ folder. Make sure the file is named 'PFXXXXX.seed' or 'PFXXXXX.full'. Delete any other files from the pfam_files/ folder. You can put any number of Pfam alignment files in the directory and sherlocc will analyse all of them.

3. Open a console or terminal window and change directory to the sherlocc directory and run sherlocc.pl.

4. You will be asked to stage at stage 3. Type 'n' and press enter.

5. Sherlocc will ask you to enter 2 parameters:

	a) A rare codon frequency threshold (i.e.: 14,15,16,etc): the lower the threshold, the more rare codon cluster regions identified will be occupied by codons with a lower frequency of usage.

	b) A window size to predict unstructured regions. We analysed all Pfam families for unstructured regions with a window size of 51. See Prilusky et al. 2005.

6. Sherlocc will ask you if you wish to clean up Stage1 and Stage2 folders after the analysis. You can decide to keep these files for further analysis (ie: the Stage 2 files for further analysis starting from stage 3).

7. You will be asked if you want Sherlocc to create the HTML outputs. If you analyse a few thousand files, the HTML output files can take some disk space. If you say no Sherlocc will still create the analysis result file.

8. The name of a directory that Sherlocc will create to store the HTML output in results/. If a directory already exists with the name you entered, Sherlocc will ask you to enter a new name or delete the content of the existing folder.


-------------------------------------
@@----If you start from stage 3----@@
-------------------------------------

1. Open a console or terminal window and change directory to the sherlocc directory and run sherlocc.pl.

2. You will be asked to stage at stage 3. Upon typing 'y', the content of Stage1 and Stage2 folders will be cleared so make sure you put the Stage2 ouput files in the Stage2/ folder only after this step. Type 'y' and press enter.

3. Put the Stage2 files in Stage2/ folder.

4. The next stepss are the same as steps 5-8 from the procedure 'Start from Stage 1.



#------------------------------------------------#
#---------What is a Sherlocc HTML output---------#
#------------------------------------------------#

The Sherlocc output displays the family protein alignment. For every protein there are two rows of information. The first row displays the amino acid, it's corresponding codon and the codon usage frequency of the specific codon. The second row prints the value calculated from the Uversky equation using the appropriate window size. It is normal not to find these values at the begining and end of the sequences as the method uses a full centered window. For a window size of 51, this window can't be calculated for the first and last 25 positions.

At the bottom of the alignment there are 3 rows: average, rare codon clusters and unstructured index values.

   a) # Average row #. Shows the codon usage frequency average at every position using a centered 7 positions wide window. For data quality purposes, windows occupied by alignment gaps (tagged as 'x') by 25% or more are excluded from the average row. The positions with a codon usage frequency average equal or lower than the specified threshold will be highlighted in orange and considered as 'slow' positions.

   b) # Rare Codon Clusters row #. This row shows the position of rare codon clusters, if any. Positions occupied by rare codon clusters are highlighted in red. To detect a rare codon cluster, Sherlocc uses a centered 7 positions wide window. Any window with 4 'slow' positions or more are considered as rare codon clusters.

   c) # Unstructured index values #. This row displays the positions of predicted family consensus unstructured regions. Positions with values below 0 are predicted to be unstructured.

Under the alignment table, two tables display the position of rare codon clusters and unstructured regions.



#------------------------------------------------#
#-------------- How to cite Sherlocc ------------#
#------------------------------------------------#

A paper has recently been submitted. Citation information will be available soon.



#------------------------------------------------#
#--------------- Contact information ------------#
#------------------------------------------------#

For additional information, questions, suggestions or bug reports please contact matthieu.chartier@usherbrooke.ca



------References------

1. Uversky VN, Gillespie JR, Fink AL (2000) Why are "natively unfolded" proteins unstructured under physiologic conditions? Proteins 41: 415-427.

2. Prilusky J, Felder CE, Zeev-Ben-Mordehai T, Rydberg EH, Man O, et al. (2005) FoldIndex: a simple tool to predict whether a given protein sequence is intrinsically unfolded. Bioinformatics 21: 3435-3438.
