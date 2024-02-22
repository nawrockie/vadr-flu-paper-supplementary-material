# 00README-tables3-4.txt
# Eric Nawrocki
# February 22, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR" related to Tables 3 and 4.
#
# =====================================================================
#
# Table 3 has data on the VADR and FLAN influenza A reference
# sequences and models.
#
# Table 4 has data on the VADR and FLAN influenza B, C and D reference
# sequences and models.
#
# There are 3 files in this directory:
# 
#  1. make-tables3-4.sh: a shell script that will download the set of
#     vadr influenza models associated with the paper and then use
#     those model files to create tables 3 and 4 from the paper by
#     calling the table-models.pl perl script.
# 
#     The first step of downloading the models will only work if you
#     have 'curl' installed. Alternatively, you can use git to clone 
#     the bitbucket repo with the models. There are instructions in
#     the 'make-tables3-4.sh' file for doing this. 
#      
#  2. clean.up.sh: a shell script that will remove all files created
#     by 'make-tables3-4.sh'.
# 
#  3. table-errors.pl: a perl script called by 'make-tables2-5-6.sh'
#     that creates 'table2.tex'.
#
#  4. table-passfail-train.pl: a perl script called by
#     'make-tables2-5-6.sh' that creates table5.tex.
#
#  5. table-passfail-test.pl: a perl script called by
#     'make-tables2-5-6.sh' that creates table6.tex.
#
#  6. parse-and-compare-flan-and-vadr-output.pl:
#  7. fetch-class-and-compare-output.pl:  
#     perl scripts called by 'make-tables2-5-6.sh' that compares flan
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
# 11. 00README-tables2-5-6.txt: this file.
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
