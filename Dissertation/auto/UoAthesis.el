(TeX-add-style-hook "UoAthesis"
 (lambda ()
    (TeX-add-symbols
     '("Chap" ["argument"] 1)
     '("ChapS" 1)
     '("acknowledgements" 1)
     '("thesisuniversity" 1)
     '("thesisfaculty" 1)
     '("thesisschool" 1)
     '("thesisyear" 1)
     '("thesispdfkeywords" 1)
     '("thesisauthor" 1)
     '("thesistitle" 1)
     '("thesisabstract" 1)
     "docspacing"
     "wordcount"
     "cleardoublepage"
     "oldchap")
    (TeX-run-style-hooks
     "etoolbox"
     "color"
     "nomencl"
     "noprefix"
     "refpage"
     "epigraph"
     "placeins"
     "subfig"
     "textpos"
     "absolute"
     "graphicx"
     "multirow"
     "footmisc"
     "bottom"
     "perpage"
     "fancyhdr"
     "times"
     "setspace"
     "pdfpages"
     "indentfirst"
     "fontenc"
     "T1"
     "url"
     "cancel"
     "makeroom"
     "bm"
     "amsmath"
     "fncychap"
     "Bjarne"
     "rep12"
     "report"
     "UKenglish"
     "a4paper"
     "12pt"
     "twoside")))

