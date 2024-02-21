# 00README-training-and-testing-sets.txt
# Eric Nawrocki
# February 21, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR" related to the training and testing sequence
# sets.
#
# =====================================================================
#
# The training and testing sets were constructed for each influenza
# type (A, B, or C) by randomly selecting sequences from two sets of
# all influenza sequences of that type available on the NCBI virus website
# (https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). The first set was
# all available sequences on March 17, 2023, and the second set was
# all available sequences released between March 18, 2023 and November
# 30, 2023. 
# 
# Training sequences and the first sets of testing sequences (referred
# to as 'test1' below) were selected from the sequences released on or
# before March 17, 2023.  The second set of testing sequences
# (referred to as 'test2' below) is more recent sequences, selected
# from the set of sequences released between March 18, 2023 and
# November 30, 2023.
#
# All selected sequences satisfied the following two criteria:
#
# - at least 60nt long (because the shortest observed sequence that
#   did not fail FLAN with a 'too short' error was 60nt).
#
# - not in the 'PAT' (patent) division, this information is available
#   in the LOCUS line of the GenBank record for each sequence 
#   (https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
# 
# The sequences were further divided into GenBank and non-GenBank
# based on the database they were deposited into. This can be
# determined from the third '|' delimited token in the sequence name: 
# e.g. : 'gi|2452516296|gb|OQ596785.1|' 
#        'gi|6406386|dbj|AB027405.1|'
# 
# If that token is 'gb' or 'ref' the sequence is eligible for the
# GenBank set. If it is 'emb' or 'dbj' the sequence is eligible for the
# non-GenBank set. 
# 
# Finally, all sets were selected to be disjoint. No sequence in a
# training set is also in a testing set and vice versa. 
#
# For flu A, 10,000 random sequences were selected 
# For flu B,  1,000 random sequences were selected 
# For flu C,    500 random sequences were selected 
# If not enough qualifying sequences existed to reach these counts,
# all qualifying sequences were selected. For example, for flu B, only
# 99 qualifying non-genbank sequences existed in the test2 dataset so
# all 99 of those were selected as the set.
#
# Below is a list of all the files in this directory, which include
# lists of all candidate sequences obtained from NCBI virus, and lists
# of the randomly selected sequences in each of the training and test
# sets for influenza type.
#
#  1. fluA.20230317.all.list: list of all candidate influenza A
#     sequences for the training and test1 testing set, released on or
#     before March 17, 2023.
#
#  2. fluA.20230318.20231130.all.list: list of all candidate influenza
#     A sequences for test2 testing set, released on or after March
#     18, 2023.
#
#  3. train.fluA.gb.10000.list: list of the randomly selected 10,000
#     sequences in the flu A genbank training set.
#
#  4. train.fluA.nongb.10000.list: list of the randomly selected
#     10,000 sequences in the flu A non-genbank training set.
#
#  5. test1.fluA.gb.10000.list: list of the randomly selected 10,000
#     sequences in the flu A genbank testing set (test1) released on
#     or before March 17, 2023.
#
#  6. test1.fluA.nongb.10000.list: list of the randomly selected
#     10,000 sequences in the flu A genbank testing set (test1)
#     released on or before March 17, 2023.
#
#  7. test2.fluA.gb.10000.list: list of the randomly selected 10,000
#     sequences in the flu A genbank testing set (test1) released on
#     or after March 18, 2023.
# 
#  8. test2.fluA.nongb.2404.list: list of the randomly selected 2,404
#     sequences in the flu A genbank testing set (test1) released on
#     or after March 18, 2023.
#
#  9. fluB.20230317.all.list: list of all candidate influenza B
#     sequences for the training and test1 testing set, released on or
#     before March 17, 2023.
#
# 10. fluB.20230318.20231130.all.list: list of all candidate influenza
#     B sequences for test2 testing set, released on or after March
#     18, 2023.
#
# 11. train.fluB.gb.1000.list: list of the randomly selected 1,000
#     sequences in the flu B genbank training set.
#
# 12. train.fluB.nongb.1000.list: list of the randomly selected
#     1,000 sequences in the flu B non-genbank training set.
#
# 13. test1.fluB.gb.1000.list: list of the randomly selected 1,000
#     sequences in the flu B genbank testing set (test1) released on
#     or before March 17, 2023.
#
# 14. test1.fluB.nongb.391.list: list of the randomly selected
#     391 sequences in the flu B genbank testing set (test1)
#     released on or before March 17, 2023.
#
# 15. test2.fluB.gb.1000.list: list of the randomly selected 1,000
#     sequences in the flu B genbank testing set (test1) released on
#     or after March 18, 2023.
# 
# 16. test2.fluB.nongb.99.list: list of the randomly selected 99
#     sequences in the flu B genbank testing set (test1) released on
#     or after March 18, 2023.
#
# 17. fluC.20230317.all.list: list of all candidate influenza C
#     sequences for the training and test1 testing set, released on or
#     before March 17, 2023.
#
# 18. fluC.20230318.20231130.all.list: list of all candidate influenza
#     C sequences for test2 testing set, released on or after March
#     18, 2023.
#
# 19. train.fluC.gb.500.list: list of the randomly selected 500
#     sequences in the flu C genbank training set.
#
# 20. train.fluC.nongb.500.list: list of the randomly selected
#     500 sequences in the flu C non-genbank training set.
#
# 21. test1.fluC.gb.16.list: list of the randomly selected 16
#     sequences in the flu C genbank testing set (test1) released on
#     or before March 17, 2023.
#
# 22. test1.fluC.nongb.500.list: list of the randomly selected
#     500 sequences in the flu C genbank testing set (test1)
#     released on or before March 17, 2023.
#
# 23. test2.fluC.gb.15.list: list of the randomly selected 15
#     sequences in the flu C genbank testing set (test1) released on
#     or after March 18, 2023.
# 
# 24. test2.fluC.nongb.130.list: list of the randomly selected 130
#     sequences in the flu C genbank testing set (test1) released on
#     or after March 18, 2023.
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
