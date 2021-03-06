\ProvidesClass{UoAthesis}

%%% Class options
\newif\if@doublesided
\@doublesidedfalse
\DeclareOption{doublesided}
{\@resitpapertrue}

\ExecuteOptions{}
\ProcessOptions

\if@doublesided
\LoadClass[UKenglish,a4paper,12pt,twoside]{report}
\else
\LoadClass[UKenglish,a4paper,12pt]{report}
\fi

%%%%%%%%%%%%%% Settable variables
\makeatletter%
\newcommand{\thesisabstract}[1]{\def \@abstract {#1}}
\newcommand{\thesistitle}[1]{\def \@thesistitle {#1}}
\newcommand{\thesisauthor}[1]{\def \@thesisauthor {#1}}
\newcommand{\thesispdfkeywords}[1]{\def \@thesispdfkeywords {#1}}
\newcommand{\thesisyear}[1]{\def \@thesisyear {#1}}
\newcommand{\thesisschool}[1]{\def \@school {#1}}
\newcommand{\thesisfaculty}[1]{\def \@college {#1}}
\newcommand{\thesisuniversity}[1]{\def \@university {#1}}
\newcommand{\acknowledgements}[1]{\def \@acknowledgements {#1}}
%%% Clear the variables at global scope
\global\let\@abstract\@empty
\global\let\@school\@empty
\global\let\@college\@empty
\global\let\@university\@empty
\global\let\@thesisyear\@empty
\global\let\@thesispdfkeywords\@empty
\global\let\@thesisauthor\@empty
\global\let\@thesistitle\@empty
\makeatother

%%% Set the default values
\thesisschool{School of Engineering}
\thesisfaculty{College of Physical Sciences}
\thesisuniversity{Univ.~of Aberdeen}
\thesispdfkeywords{}

%%% Other details
\newcommand{\docspacing}{\onehalfspacing}
%%% The university guidelines state 25mm top and bottom margins, and
%%% 30mm left and right margins. No binding indent is allocated
\newlength{\bindingedgeindent}\setlength{\bindingedgeindent}{0mm}
\newlength{\sideMargin}\setlength{\sideMargin}{30mm}
\newlength{\topMargin}\setlength{\topMargin}{25mm}
\newlength{\bottomMargin}\setlength{\bottomMargin}{25mm}
\RequirePackage[a4paper,bindingoffset=\bindingedgeindent,%
hdivide={\sideMargin,*,\sideMargin},%
vdivide={\topMargin,*,\bottomMargin}]{geometry}

%%% Chapter Formatting
\RequirePackage[Bjarne]{fncychap} %% Chapter style

%%% Reference sections
\RequirePackage[backend=bibtex,style=chem-acs, doi=true, articletitle=true, %
biblabel=brackets,refsection=chapter]{biblatex}
%\RequirePackage[style=nature, article-title=true]{biblatex}
\defbibheading{thesisChapterBib}{\section{References}}

%%% Standard latex packages
\RequirePackage{amsmath}
\RequirePackage{bm}
\RequirePackage[makeroom]{cancel}
\RequirePackage{url}
\RequirePackage[T1]{fontenc} %Better fonts
\RequirePackage{indentfirst} %indent first line of sections
\RequirePackage{pdfpages} %% Allows insertion of pdfs as pages
\RequirePackage{setspace} % Document spacing
\RequirePackage{times}

%%% Header/Footer Format
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[RO]{\slshape \rightmark}
\fancyhead[LE]{\slshape \leftmark}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\addtolength{\headheight}{2.5pt}
\RequirePackage[bottom,perpage]{footmisc} %% Footnotes at bottom of page

%%% Tables
\RequirePackage{multirow} %% Allows multiple row items in tables

%%% Make equations numbered by section
\numberwithin{equation}{section}

%%% Image manipulation
\usepackage{graphicx}
\usepackage[absolute]{textpos} %for the cover image
\setlength{\TPHorizModule}{1mm} %sets the textpos unit
\setlength{\TPVertModule}{\TPHorizModule} 
\usepackage{subfig} %for multi image file floats
%%% This command lets figures fill 85 percent of the top of a
%%% page without placing it on its own
\renewcommand{\topfraction}{0.85}
\renewcommand{\textfraction}{0.1}
%%% Makes sure a page with 75
\renewcommand{\floatpagefraction}{0.75}

%%% This makes clear double page not print headers or footers
\makeatletter
\def\cleardoublepage{\clearpage\if@twoside \ifodd\c@page\else
\hbox{}
\vspace{\fill}
\thispagestyle{plain}
\newpage
\if@twocolumn\hbox{}\newpage\fi\fi\fi}
\makeatother

%%% Keep floats in the chapter in which they were placed
\usepackage{placeins}
\let\oldchap=\chapter
\renewcommand*{\chapter}{%
  \secdef{\Chap}{\ChapS}% whether there's a star or not...
}
\newcommand\ChapS[1]{\FloatBarrier\singlespacing\cleardoublepage\oldchap*{{\uppercase{#1}}}\markboth{{\uppercase{#1}}}{{\uppercase{#1}}}\docspacing}
\newcommand\Chap[2][]{\FloatBarrier\singlespacing\cleardoublepage\oldchap[#1]{{\uppercase{#2}}}\docspacing}

%% Epigraphs for stating where chapters were from
\RequirePackage{epigraph}
\setlength{\epigraphwidth}{0.9\columnwidth}
\setlength{\epigraphrule}{0pt}
\setlength{\beforeepigraphskip}{-\baselineskip}
\setlength{\afterepigraphskip}{\baselineskip}
\renewcommand{\textflush}{flushright}
\renewcommand{\epigraphsize}{\footnotesize}

%%% Nomenclature
\RequirePackage{ifthen}
\RequirePackage[noprefix,refpage]{nomencl} %My version of the Nomencl package
\setlength{\nomitemsep}{-\parsep} %single line spacing nomencl
\renewcommand{\nomgroup}[1]{%
 \ifthenelse{\equal{#1}{A}}{\vspace{10pt} \item[\textbf{Acronyms/Terminology}]}{%
 \ifthenelse{\equal{#1}{F}}{\vspace{10pt} \item[\textbf{Functions}]}{%
 \ifthenelse{\equal{#1}{O}}{\vspace{10pt} \item[\textbf{Operators}]}{%
 \ifthenelse{\equal{#1}{S}}{\vspace{10pt} \item[\textbf{Notation}]}{%
 \ifthenelse{\equal{#1}{V}}{\vspace{10pt} \item[\textbf{Variables/Constants}]}{}}}}}}

%%% Stop links being colored in the pdf, place
\RequirePackage{color}
\definecolor{Navy}{rgb}{0.1,0.1,0.45}
\makeatletter
\RequirePackage[pdftex,%
pdftitle={\@thesistitle},
plainpages=false,
pdfpagelabels,
hypertexnames=true,
naturalnames,
pdfauthor={\@thesisauthor},
pdfsubject={PhD Thesis},
pdfkeywords={\@thesispdfkeywords},
pdfproducer={Latex},
pdfcreator={pdflatex},
bookmarks,bookmarksnumbered,
%pagebackref,
colorlinks,
linkcolor=black,
citecolor=black,
filecolor=black,
urlcolor=black,
breaklinks=true
]{hyperref}
\makeatother

%%% Allow hooking into the start of the document
\RequirePackage{etoolbox}

%%% Makes page numbers hyperlinked
\renewcommand*{\pagedeclaration}[1]{\unskip, \hyperpage{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Fix the colours of transparent images
\pdfpageattr {/Group << /S /Transparency /I true /CS /DeviceRGB>>}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Make it look like a book in acroread
\makeatletter
\if@doublesided
\pdfcatalog{/PageLayout /TwoPageRight}
\fi
\makeatother

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Command to collect nomenclature data
\makenomenclature
\makeatletter
\AfterEndPreamble{%
  \singlespacing
  %%% Title page
  \ifpdf\phantomsection\pdfbookmark[0]{Front Cover}{FrontCover}\fi
  \begin{titlepage}
    \begin{textblock}{63}(25,30)
      \includegraphics[clip,width=\textwidth]{figures/UoALogo}
    \end{textblock}
    \mbox{}\vskip50mm%
    \begin{center}
      \Large\expandafter\textsc\expandafter{\@thesistitle}
    \end{center}
    \vfill
    \begin{center}
      \rm \textbf{\textsc\@thesisauthor}\\
      \@thesisyear
    \end{center}
    \vskip1.5in%
    \begin{center}
      \large \textsc{ This thesis was submitted as part of the
        requirement for the MEng.~degree in Engineering.}
    \end{center}
    \vfill
    \begin{center}
      \@school\hfill\@thesisyear\hfill\@university
    \end{center}
  \end{titlepage}

  %%% Abstract
  \pagenumbering{roman}
  \setcounter{page}{1}
  \chapter*{Abstract}
  \addcontentsline{toc}{chapter}{Abstract}%
  \@abstract
  \cleardoublepage

  %%% Table of Contents
  \cleardoublepage
  \ifpdf\phantomsection\pdfbookmark[0]{Contents}{Contents}\fi
  \tableofcontents
  \vfill
  %%% Required at the end of the table of contents
  \newcommand{\wordcount}{RUN MakeFigsBibIndex}
  \InputIfFileExists{wc}{}{}
  \begin{center}
    {\bf Final word count:} \wordcount
  \end{center}
  \cleardoublepage
  \ifpdf\phantomsection\pdfbookmark[0]{List of Figures}{TableoffiguresPg}\fi
  \listoffigures
  \cleardoublepage
  \ifpdf\phantomsection\pdfbookmark[0]{List of Tables}{TableoftablesPg}\fi
  \listoftables
  \cleardoublepage
  \phantomsection\addcontentsline{toc}{chapter}{Nomenclature}
  \printnomenclature
  \chapter*{Acknowledgements}
  \addcontentsline{toc}{chapter}{Acknowledgements}
  \@acknowledgements
  \cleardoublepage
  \docspacing
  \pagenumbering{arabic}
}
\makeatother
