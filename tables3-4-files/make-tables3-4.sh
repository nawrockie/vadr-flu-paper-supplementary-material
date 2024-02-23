#!/bin/bash
VERSION=1.6.3-2dev1

if [ ! -d "../vadr-models-flu-$VERSION" ]; then
    if [ -e "../vadr-models-flu-$VERSION.tar.gz" ]; then
        echo "ERROR: ../vadr-models-flu-$VERSION/ directory does not exist."
        echo "Execute the following 'cd ..; tar xf vadr-models-flu-$VERSION.tar.gz;'"
        echo "Then return to this directory and try again."
    else 
        echo "ERROR ../vadr-models-flu-$VERSION/ directory does not exist."
        echo "and ../vadr-models-flu-$VERSION.tar.gz does not exist."
        echo "Try redownloading the supplementary data and trying again."
    fi
    exit 1;
fi

# create the table
perl table-models.pl ../vadr-models-flu-$VERSION/flu.minfo ../vadr-models-flu-$VERSION/flan ../vadr-models-flu-$VERSION > tables3-4.tex

