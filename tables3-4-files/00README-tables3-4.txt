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
#  3. table-models.pl: a perl script called by 'make-tables3-4.sh' 
#     that creates the file 'tables3-4.tex'
# 
#  4. 00README-tables3-4.txt: this file.
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
