# 00README-fig1-nseq.txt
# Eric Nawrocki
# June 13, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR" related to Figure 1
#
# =====================================================================
#
# Figure 1 shows how the number of influenza sequences in INSDC
# databases has changed between 2000 and 2024. There are 4 files in
# this directory, described below:
#
# 1. nseq.data: the raw counts obtained from the NCBI Virus resource
#    on June 13, 2024 by restricting 'release date' to between
#    1/1/1900 and 12/31 of each year. The following taxids were used
#    for each influenza type: A: 11320, B: 11520, C: 11552, D:
#    1511084. 
#
# 2. fig-nseq-data2R.pl: perl script that takes nseq.data as input and
#    creates an R script that will create a PDF similar to figure 1.
# 
# 3. make-nseq-R-script.sh: shell script that will call
#    fig-nseq-data2R.pl to make the R script and then execute the R
#    script to make a PDF.
#
# 4. 00README-fig1-nseq.txt: this file
#
################################################
# Question/problems? email eric.nawrocki@nih.gov

