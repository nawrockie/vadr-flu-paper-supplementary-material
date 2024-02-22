#!/bin/bash
VERSION=1.6.3-2dev1

############################################
curl -k -L -o vadr-models-flu-$VERSION.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu/$VERSION/vadr-models-flu-$VERSION.tar.gz 
tar xf vadr-models-flu-$VERSION.tar.gz 
############################################

############################################
# if above command doesn't work, comment out that line and uncomment the block below
# to get the models by cloning the bitbucket repository for the models and checking 
# out the version associated with the paper (1.6.3-2dev1)
#
####beginning of block####
#git clone https://bitbucket.org/nawrockie/vadr-models-flu
#mv vadr-models-flu vadr-models-flu-$VERSION
#cd vadr-models-flu-$VERSION
#git checkout vadr-models-flu-$VERSION
#rm -rf .git
#cd ..
####end of block####
#########################################

# create the table
perl table-models.pl vadr-models-flu-$VERSION/flu.minfo vadr-models-flu-$VERSION/flan vadr-models-flu-$VERSION > tables3-4.tex

