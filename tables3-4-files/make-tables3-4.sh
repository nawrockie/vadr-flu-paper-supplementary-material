#!/bin/bash
VERSION=1.6.3-2dev1

cd ..
tar xf vadr-models-flu-$VERSION.tar.gz
cd tables3-4-files

# create the table
perl table-models.pl vadr-models-flu-$VERSION/flu.minfo ../vadr-models-flu-$VERSION/flan ../vadr-models-flu-$VERSION > tables3-4.tex

