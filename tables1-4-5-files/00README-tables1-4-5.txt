# 00README-tables2-5-6.txt
# Eric Nawrocki
# February 20, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR" related to Tables 1, 4 and 5.
#
# =====================================================================
#
# Table 1 has data on how consistent FLAN and VADR errors are in the
# combined training and testing sets.
#
# Table 4 has data on a comparison of the pass/fail outcomes of the
# training set sequences from FLAN and VADR.
#
# Table 5 has data on a comparison of the pass/fail outcomes of the
# testing set sequences from FLAN and VADR.
#
# The 11 files in this directory are necessary to reproduce tables 1, 4,
# and 5 using data in the top-level 'flan-train-test-output' and
# 'vadr-train-test-output' directories one level up from this
# directory.
# 
# A brief description of each of the 11 files is below:
# 
#  1. make-tables1-4-5.sh: a shell script that will create more than
#     200 files, most importantly table1.tex, table4.tex and
#     table5.tex which are the same tables in the paper with the
#     exception that the table 1 in the paper is a manually edited
#     version of table1.tex with formatting changes and reordering of
#     the rows. This script will also create the two files
#     'sum.coords.txt' and 'sum.info.txt' which include data on how
#     many of the predictions between vadr and flan were identical or
#     different, and how many of the classifcations (type and subytpe)
#     were identical or different, respectively. 
#     Execute this file with the command 'sh ./make-tables1-4-5.sh'
#
#  2. clean.up.sh: a shell script that will remove all files created
#     by 'make-tables1-4-5.sh'.
# 
#  3. table-errors.pl: a perl script called by 'make-tables1-4-5.sh'
#     that creates 'table1.tex'.
#
#  4. table-passfail-train.pl: a perl script called by
#     'make-tables1-4-5.sh' that creates table4.tex.
#
#  5. table-passfail-test.pl: a perl script called by
#     'make-tables1-4-5.sh' that creates table5.tex.
#
#  6. parse-and-compare-flan-and-vadr-output.pl:
#  7. fetch-class-and-compare-output.pl:  
#     perl scripts called by 'make-tables1-4-5.sh' that compares flan
#     and vadr output and creates intermediate files used by the other
#     perl scripts
# 
#  8. coords.pl: a perl script that summarizes the per feature
#     coordinate comparisons between flan and vadr.
#
#  9. info.pl: a perl script that summarizes the comparison of the
#     classification (type/subtype) of each sequence between flan and
#     vadr.
#  
# 10. product2gene.map.txt: a map of product names to gene names, used
#     by 'parse-and-compare-flan-and-vadr-output.pl' when parsing flan
#     output. 
#
# 11. 00README-tables1-4-5.txt: this file.
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
