# create the R script 
perl fig-nseq-data2R.pl nseq.data > nseq.R
# execute the R script (requires R to be in your path)
R --vanilla < nseq.R
# will create nseq.pdf
