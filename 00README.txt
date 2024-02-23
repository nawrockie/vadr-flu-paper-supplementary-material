# 00README.txt
# Eric Nawrocki
# February 20, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR".
#
# Files that end in .gz can be unpacked with the command 'gunzip
# filename'
# 
# =====================================================================
#
# There are six subdirectories in this directory, each contains its
# own 00README file that explains the files in that subdirectory. In
# each of the three subdirectories with names beginning with 'table'
# there are shell scripts that can be executed to reproduce the tables
# in the paper. Each subdirectory is described below:
#
# train-and-test-sets/:   lists of sequences in the training and testing
#                         sequence sets
#
# flan-train-test-output: FLAN output files for the training and
#                         testing sets
#
# vadr-train-test-output: VADR output files for the training and
#                         testing sets
# 
#
# table1-files/:          files relevant to table 1 which has counts
#                         of influenza sequences deposited in GenBank
#                         since 2018.
#
# tables2-5-6-files/:     files relevant to tables 2, 5, and 6 which have
#                         pass/fail outcome and error comparisons
#                         between FLAN and VADR on the training and
#                         testing datasets.
#
# tables3-4-files/:       files relevant to tables 3 and 4 which have
#                         information on reference and model sequences
#                         used by FLAN and VADR.
#
# There are two additional files:
# 
# vadr-models-flu-1.6.3-2.tar.gz: G'zipped tarball of the vadr 1.6.3-2
#                         models used by VADR for the paper. To unpack
#                         this execute:
#                         'tar xf vadr-models-flu-1.6.3-2.tar.gz'
#                         You will need to unpack this file before you 
#                         can run the 'make-tables3-4.sh' script in 
#                         'tables3-4-files/".
#                         Unpacking will create a directory named: 
#                         'vadr-models-flu-1.6.3-2' which will include 
#                         a '00NOTES.txt' and '00README.txt' files
#                         with more information on the models as well 
#                         as a 'mapping-flan-to-genbank/' directory
#                         with information on how the FLAN reference
#                         sequences were mapped to GenBank sequences.
# 
# 00README.txt:           this file
#
# ====================================================================
# Question/problems/feature requests? email eric.nawrocki@nih.gov


