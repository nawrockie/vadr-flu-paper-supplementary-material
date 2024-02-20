# 00README-table1.txt
# Eric Nawrocki
# February 20, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR" related to Table 1.
#
# =====================================================================
#
# Table 1 shows the number of influenza sequences published in GenBank
# with release dates in each year between 2018 and 2023.
# 
# This directory contains 32 list files with names ending in either
# '.gblist' or '.list': 
#
#  - 24 of the '.gblist' files contain the accessions published in GenBank for each of the six
#    years and four types of influenza (A, B, C, and D). An example is 
#    'fluC.2022.gblist' which contains all fluC sequences published in
#    GenBank in 2022. 
# 
#  - 4 of the remaining 8 files end in '.gblist' and include all GenBank
#    accessions (from any date) for each type as of February 20, 2024.
#    An example is 'fluB.all.20240220.gblist' which includes all flu B
#    GenBank accessions. 
# 
#  - The final 4 files end in '.list' and include all INSDC (GenBank,
#    ENA and DDBJ) accessions from any date as of February 20,
#    2024. An example is 'fluA.all.20240220.list' which includes all
#    flu A accesions.
#
# These lists were obtained on February 20, 2024 from the NCBI Virus
# interactive dashboard, tabular view: 
# https://www.ncbi.nlm.nih.gov/labs/virus/vssi
# by searching for either "influenza A virus" (taxid 11320),
# "influenza B virus" (taxid 11520), "influenza C virus" (taxid 11552)
# or "influenza D virus (taxid 1511084) and restricting by 'release
# date'. 
#
# The dashboard allows you to download all INSDC accessions, which
# include sequences originally published in GenBank, ENA and DDBJ.
# For the 6 per-year rows, these were filtered to only those
# accessions published in GenBank by fetching all the sequences using
# an internal NCBI script that fetches sequences and renames them with
# a 'long form' sequence name that includes the GI and the original
# database, e.g.
#
# gi|117572924|gb|CY017309.1|
#
# These sequence files were filtered for sequence names containing
# 'gb' to restrict to only sequences deposited in GenBank.
# 
# The lists in this directory include only the resulting GenBank
# accessions. 
#
# The accession list from NCBI virus can also be split into GenBank vs
# ENA/DDBJ based on the first two letters of the accessions. The
# 481,079 GenBank 2018 to 2023 influenza sequences from Table 1 all
# begin with one of the following two letter prefixes: 
# CY, JQ, KJ, KM, KP, KR, KT, KU, KX, KY, MF, MG, MH, MK, MN, MT, MW,
# MZ, OK, OL, OM, ON, OP, OQ, OR, PP.
# 
# See https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
# for more information on accession prefixes.
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
