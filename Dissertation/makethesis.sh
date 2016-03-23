#/bin/bash
cd figures; ./Makefigures; cd ..
pdflatex Dissertation.tex --draftmode || exit 1

makeindex Dissertation.nlo -s nomencl.ist -o Dissertation.nls  || exit 1
bibtex Dissertation  || exit 1
for file in $(ls Dissertation*-blx.aux); do
    bibtex $file 
done
Words=$(pdftotext Dissertation.pdf - | wc -w)

WC=$(./texcount/texcount.pl -inc -1 -sum=1,1,1 -q Dissertation.tex | gawk '{print $1}')

pdflatex Dissertation.tex --draftmode || exit 1
pdflatex Dissertation.tex  || exit 1

echo "*******************************************"
echo "       Estimated Word Count = $WC"
echo "*******************************************"

echo "\renewcommand{\wordcount}{$WC}" > wc.tex