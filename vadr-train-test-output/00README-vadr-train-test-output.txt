# 00README-vadr-train-test-output.txt
# Eric Nawrocki
# February 22, 2024
# 
# This file explains the contents of this directory which contains
# supplementary data for the paper "Influenza sequence validation and
# annotation using VADR", specifically output files from VADR 1.6.3 
# run on the training and testing sets described in the paper.
#
# =====================================================================
#
# There are 156 files in this directory, these are all output from VADR
# 1.6.3 when run on the training and testing sets, but not all VADR
# output files are included here. There are six types of output files
# here, which are used by the scripts in the 'tables2-5-6-files'
# directory to create tables 2, 5 and 6 in the paper.
#
# Six types of vadr output files:
# .vadr.pass.list: list of all passing sequences
# .vadr.fail.list: list of all failing sequences
# .vadr.pass.tbl:  feature table for all passing sequences
# .vadr.fail.tbl:  feature table for all failing sequences
# .vadr.sqa:       tabular per sequence file with information on the
#                  annotation of each sequence
# .vadr.alt:       tabular per alert file with information on each
#                  alert reported by vadr
#
# More information on the format of these files can be found here:
# https://github.com/ncbi/vadr/blob/vadr-1.6.3/documentation/formats.md
#
# These six sets of files exist for 18 different runs of VADR
# (v-annotate.pl v1.6.3) using recommended options and the vadr flu
# models version 'vadr-models-flu-1.6.3-2' with the command: 
#
# v-annotate.pl -f -r --mdir $VADRFLUDIR --mkey flu --atgonly \
# --alt_fail extrant5,extrant3 --xnocomp --nomisc <fastafile> <outdir>
# 
# where $VADRFLUDIR is the path to the vadr-models-flu-1.6.3-2
# directory. These are the 4*18=72 files that do not include '.fo' or
# '.nto' in their names. An example of one of these 72 files is:
# train.fluA.gb.10000.vadr.sqa
# 
# There are 12 additional sets of 4 files, created by using different
# models (the 'FLAN' and 'FLAN-ntonly' rows in table 5). The 6 sets of
# files with '.fo.' in their names correspond to the 'FLAN' rows and
# were created with this command:
# 
# v-annotate.pl -f -r --mdir $VADRFLUDIR/flan --mkey flu-flan --atgonly \
# --alt_fail extrant5,extrant3 --xnocomp --nomisc <fastafile> <outdir>
# 
# The 6 sets of 4 files with '.nto.' in their names correspond with
# 'FLAN-ntonly' and were created with this
# command: 
#
# v-annotate.pl -f -r --mdir $VADRFLUDIR/flan-nt -i \
# $VADRFLUDIR/flan-nt/flu-flan-nt.minfo -m \
# $VADRFLUDIR/flan/flu-flan.cm -n $VADRFLUDIR/flan/flu-flan.fa -x \
# $VADRFLUDIR/flan-nt --atgonly \
# --alt_fail extrant5,extrant3 --xnocomp --nomisc <fastafile> <outdir>
#
################################################
# Question/problems? email eric.nawrocki@nih.gov
