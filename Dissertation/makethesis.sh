#/bin/bash
pdflatex main.tex --draftmode || exit 1

makeindex main.nlo -s nomencl.ist -o main.nls  || exit 1
bibtex main  || exit 1
for file in $(ls main*-blx.aux); do
    bibtex $file 
done
Words=$(pdftotext main.pdf - | wc -w)

WC=$(./texcount/texcount.pl -inc -1 -sum=1,1,1 -q main.tex | gawk '{print $1}')

pdflatex main.tex --draftmode || exit 1
pdflatex main.tex  || exit 1

echo "*******************************************"
echo "       Estimated Word Count = $WC"
echo "*******************************************"

echo "\renewcommand{\wordcount}{$WC}" > wc.tex