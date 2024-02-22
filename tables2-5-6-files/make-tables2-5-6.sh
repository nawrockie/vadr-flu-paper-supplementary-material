set -e
# remove old files
touch sum.passfail.txt
touch sum.info.txt
touch sum.coords.txt
rm sum.passfail.txt
rm sum.info.txt
rm sum.coords.txt

# make table2.data.txt
for a in \
train.fluA.gb.10000    \
train.fluA.nongb.10000 \
train.fluB.gb.1000     \
train.fluB.nongb.1000  \
train.fluC.gb.500      \
train.fluC.nongb.500   \
test1.fluA.gb.10000    \
test1.fluA.nongb.10000 \
test1.fluB.gb.1000     \
test1.fluB.nongb.391   \
test1.fluC.gb.16       \
test1.fluC.nongb.500   \
test2.fluA.gb.10000    \
test2.fluA.nongb.2404  \
test2.fluB.gb.1000     \
test2.fluB.nongb.99    \
test2.fluC.gb.15       \
test2.fluC.nongb.130   \
; do
    perl parse-and-compare-flan-and-vadr-output.pl --deput --misc2cds \
        flan-output/$a.flan.tbl \
        vadr-output/$a.vadr.pass.tbl \
        vadr-output/$a.vadr.fail.tbl \
        vadr-output/$a.vadr.alt \
        vadr-output/$a.vadr.sqa \
        product2gene.map.txt \
        compare.$a > \
        compare.$a.txt
    for c in FPVP FPVF FFVP FFVF; do 
        perl fetch-class-from-compare-output.pl compare.$a.txt $c > compare.$c.$a.txt
        # write to passfail summary file
        echo -n $a " " $c " " >> sum.passfail.txt
        grep SEQUENCE compare.$c.$a.txt | wc -l >> sum.passfail.txt
    done
    # write to info summary file
    for c in identical different; do 
        echo -n $a " " $c " " >> sum.info.txt
        grep SEQUENCE compare.$a.txt | grep info:$c | wc -l >> sum.info.txt
        echo "" >> sum.info.txt
    done
    # write to coords summary file
    for c in coords-same coords-diff vadronly flanonly; do
        echo -n $a " " $c " " >> sum.coords.txt
        grep FEATURE compare.$a.txt | grep $c | grep -v "\:$c" | wc -l >> sum.coords.txt
        echo "" >> sum.coords.txt
    done
done
# combine all flan.err.list and vadr.fatal.list files made by parse-and-compare-flan-and-vadr-output.pl above:
cat compare.train.*.flan.err.list   compare.test1.*.flan.err.list   compare.test2.*.flan.err.list   > all.flan.err.list
cat compare.train.*.vadr.fatal.list compare.test1.*.vadr.fatal.list compare.test2.*.vadr.fatal.list > all.vadr.fatal.list
# make table2.tex
perl ./table-errors.pl parse-and-compare-flan-and-vadr-output.pl all.flan.err.list all.vadr.fatal.list > table2.tex
# table 2 in the paper was manually edited from the 'table2.tex' file created above to fix spacing and reorder some errors

# make table5.data.txt
touch table5.data.txt
rm table5.data.txt
for a in \
train.fluA.gb.10000    \
train.fluA.nongb.10000 \
train.fluB.gb.1000     \
train.fluB.nongb.1000  \
train.fluC.gb.500      \
train.fluC.nongb.500   \
; do
    perl parse-and-compare-flan-and-vadr-output.pl --deput --misc2cds \
        flan-output/$a.flan.tbl \
        vadr-output/$a.fo.vadr.pass.tbl \
        vadr-output/$a.fo.vadr.fail.tbl \
        vadr-output/$a.fo.vadr.alt \
        vadr-output/$a.fo.vadr.sqa \
        product2gene.map.txt \
        fo.compare.$a > \
        fo.compare.$a.txt
    perl parse-and-compare-flan-and-vadr-output.pl --deput --misc2cds \
        flan-output/$a.flan.tbl \
        vadr-output/$a.nto.vadr.pass.tbl \
        vadr-output/$a.nto.vadr.fail.tbl \
        vadr-output/$a.nto.vadr.alt \
        vadr-output/$a.nto.vadr.sqa \
        product2gene.map.txt \
        nto.compare.$a > \
        nto.compare.$a.txt
    for c in FPVP FPVF FFVP FFVF; do 
        perl fetch-class-from-compare-output.pl fo.compare.$a.txt  $c > fo.compare.$c.$a.txt
        perl fetch-class-from-compare-output.pl nto.compare.$a.txt $c > nto.compare.$c.$a.txt
        # write to passfail summary file
        echo -n $a " " $c " " fin " " >> table5.data.txt
        # the following file was created in the first block
        grep SEQUENCE compare.$c.$a.txt | wc -l >> table5.data.txt
        echo -n $a " " $c " " fo " " >> table5.data.txt
        grep SEQUENCE fo.compare.$c.$a.txt | wc -l >> table5.data.txt
        echo -n $a " " $c " " nto " " >> table5.data.txt
        grep SEQUENCE nto.compare.$c.$a.txt | wc -l >> table5.data.txt
    done
done
perl table-passfail-train.pl table5.data.txt > table5.tex

# make table6.data.txt
touch table6.data.txt
rm table6.data.txt
# make table2.data.txt
for a in \
test1.fluA.gb.10000    \
test1.fluA.nongb.10000 \
test1.fluB.gb.1000     \
test1.fluB.nongb.391   \
test1.fluC.gb.16       \
test1.fluC.nongb.500   \
test2.fluA.gb.10000    \
test2.fluA.nongb.2404  \
test2.fluB.gb.1000     \
test2.fluB.nongb.99    \
test2.fluC.gb.15       \
test2.fluC.nongb.130   \
; do
    for c in FPVP FPVF FFVP FFVF; do 
        echo -n $a " " $c " " fin " " >> table6.data.txt
        grep SEQUENCE compare.$c.$a.txt | wc -l >> table6.data.txt
    done
done
perl table-passfail-test.pl table6.data.txt > table6.tex
