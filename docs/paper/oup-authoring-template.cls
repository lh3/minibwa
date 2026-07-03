%%
%% Simplified OUP authoring template for [webpdf,contemporary,large,namedate]
%% Based on OUP's original oup-authoring-template.cls Version 1.1, 2022
%%

\newcommand\classname{oup-authoring-template}
\newcommand\lastmodifieddate{2022/03/18}
\newcommand\versionnumber{1.1}

\NeedsTeXFormat{LaTeX2e}[2001/06/01]
\ProvidesClass{\classname}[\lastmodifieddate\space\versionnumber]

\newif\if@cropmarkson \@cropmarksontrue
\newif\if@numbib \@numbibtrue
\newif\if@unnumsec\global\@unnumsecfalse
\newif\if@final\global\@finalfalse
\newif\if@weblink\@weblinktrue
\newif\if@proof\global\@prooffalse

\if@compatibility\else
\DeclareOption{namedate}{\PassOptionsToPackage{authoryear,round}{natbib}\global\@numbibfalse}
\DeclareOption{draft}{\PassOptionsToPackage{draft}{graphicx}}
\DeclareOption{b4paper}{\PassOptionsToPackage{b4}{crop}}
\DeclareOption{centre}{\PassOptionsToPackage{center}{crop}}
\DeclareOption{crop}{\PassOptionsToPackage{cam}{crop}\global\@cropmarksontrue}
\DeclareOption{nocrop}{\PassOptionsToPackage{off}{crop}\global\@cropmarksonfalse}
\DeclareOption{info}{\PassOptionsToPackage{info}{crop}}
\DeclareOption{noinfo}{\PassOptionsToPackage{noinfo}{crop}}
\DeclareOption{final}{\global\@finaltrue}
\DeclareOption{unnumsec}{\global\@unnumsectrue}
\DeclareOption{webpdf}{\@weblinktrue}
\DeclareOption{modern}{}
\DeclareOption{traditional}{}
\DeclareOption{contemporary}{}
\DeclareOption{large}{}
\DeclareOption{medium}{}
\DeclareOption{small}{}
\fi
\ExecuteOptions{b4paper,centre,info}
\ProcessOptions

% Contemporary large paper size
\setlength{\paperheight}{275.84truemm}
\setlength{\paperwidth}{209.80truemm}

% Load all necessary packages
\RequirePackage{crop}
\RequirePackage{graphicx}
\RequirePackage{caption}
\RequirePackage{amsmath}
\RequirePackage{array}
\RequirePackage{color}
\RequirePackage{xcolor}
\RequirePackage{amssymb}
\RequirePackage{flushend}
\RequirePackage{stfloats}
\RequirePackage[figuresright]{rotating}
\RequirePackage{chngpage}
\RequirePackage{totcount}
\RequirePackage{fix-cm}

\def\sffamilyfont{\sffamily}
\def\sffamilyfontitalic{\sffamily\itshape\selectfont}
\def\sffamilyfontbold{\sffamily\bfseries\selectfont}
\def\sffamilyfontbolditalic{\sffamily\bfseries\itshape\selectfont}
\def\sffamilyfontcn{\sffamily\fontseries{m}\fontshape{n}\selectfont}
\def\sffamilyfontcnitalic{\sffamily\fontseries{m}\fontshape{it}\selectfont}
\def\sffamilyfontcnbold{\sffamily\bfseries\selectfont}
\def\sffamilyfontcnbolditalic{\sffamily\fontseries{b}\fontshape{it}\selectfont}

\newcommand\@ptsize{0}
\@twosidetrue
\@mparswitchfalse
\reversemarginpar

% Contemporary large normal size
\renewcommand\normalsize{%
   \@setfontsize\normalsize{7.5bp}{10.5pt}%
   \abovedisplayskip 11\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\normalsize
\let\@bls\baselineskip

\newcommand\small{%
    \@setfontsize\small{7}{10}%
    \abovedisplayskip 10\p@ minus 3\p@
    \belowdisplayskip \abovedisplayskip
    \abovedisplayshortskip \z@ plus 2\p@
    \belowdisplayshortskip 4\p@ plus 2\p@ minus2\p@
    \def\@listi{\topsep 4.5\p@ plus 2\p@ minus 1\p@
       \itemsep \parsep
       \topsep 4\p@ plus 2\p@ minus 2\p@}}

\newcommand\footnotesize{%
    \@setfontsize\footnotesize{8}{10}%
    \abovedisplayskip 6\p@ minus 3\p@
    \belowdisplayskip\abovedisplayskip
    \abovedisplayshortskip \z@ plus 3\p@
    \belowdisplayshortskip 6\p@ plus 3\p@ minus 3\p@
    \def\@listi{\topsep 3\p@ plus 1\p@ minus 1\p@
       \parsep 2\p@ plus 1\p@ minus 1\p@\itemsep \parsep}}

\def\scriptsize{\@setfontsize\scriptsize{6.5pt}{9.5pt}}
\def\tiny{\@setfontsize\tiny{5pt}{7pt}}
\def\large{\@setfontsize\large{11.5pt}{12pt}}
\def\Large{\@setfontsize\Large{14pt}{16}}
\def\LARGE{\@setfontsize\LARGE{15pt}{17pt}}
\def\huge{\@setfontsize\huge{22pt}{22pt}}
\def\Huge{\@setfontsize\Huge{30pt}{30pt}}

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamilyfont}{\mathsf}
\DeclareOldFontCommand{\sfit}{\normalfont\sffamily\itshape}{\mathsf}
\DeclareOldFontCommand{\sfb}{\normalfont\sffamilyfontbold}{\mathsf}
\DeclareOldFontCommand{\sfbi}{\normalfont\sffamily\bfseries\itshape}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}

\def\oddsideskip{48pt}
\def\evensideskip{56pt}

\newdimen\croppaperwidth
\newdimen\croppaperheight
\setlength{\croppaperwidth}{\paperwidth}
\setlength{\croppaperheight}{\paperheight}
\CROP@size{\croppaperwidth}{\croppaperheight}

\newdimen\Croppdfwidth
\newdimen\Croppdfheight
\newdimen\Trimpdfwidth
\newdimen\Trimpdfheight
\Croppdfwidth=\croppaperwidth
\Croppdfheight=\croppaperheight
\advance\Croppdfwidth by -0.71mm
\advance\Croppdfheight by -0.92mm
\Trimpdfwidth=\paperwidth
\Trimpdfheight=\paperheight
\advance\Trimpdfwidth by -0.59mm
\advance\Trimpdfheight by -0.88mm

\newdimen\CP@toff@wd
\newdimen\CP@toff@ht
\newdimen\CP@boff@wd
\newdimen\CP@boff@ht
\newdimen\CP@crop@wd
\newdimen\CP@crop@ht
\newdimen\CP@bled@wd
\newdimen\CP@bled@ht
\newdimen\CP@trim@wd
\newdimen\CP@trim@ht

\def\str@yes{yes}
\def\SetCrop#1#2{%
    \gdef\IsCropSet{yes}
    \global\CP@crop@wd=#1\relax
    \global\CP@crop@ht=#2\relax}
\def\SetTrim#1#2{%
    \gdef\IsTrimSet{yes}
    \global\CP@trim@wd=#1\relax
    \global\CP@trim@ht=#2\relax}
\def\SetBleed#1#2{%
    \gdef\IsBleedSet{yes}
    \global\CP@bled@wd=#1\relax
    \global\CP@bled@ht=#2\relax}

\everyjob\expandafter{%
    \the\everyjob
    \typeout{%
        \filename\space <\filedate>^^J
        Version: v\fileversion^^J
        LaTeX macros for setting Page Box parameters
    }%
  \IfFileExists{\filename.cfg}{%
        \begingroup\@@input\@filef@und\endgroup
    }{%
        \typeout{No File: \filename.cfg}%
    }%
  \IfFileExists{\jobname.cfg}{%
        \begingroup\@@input\@filef@und\endgroup
    }{%
        \typeout{No File: \jobname.cfg}%
    }%
}

\SetCrop{\Trimpdfwidth}{\Trimpdfheight}
\SetTrim{\Trimpdfwidth}{\Trimpdfheight}
\SetBleed{0mm}{0mm}

\def\do@pagebox@calc{%
    \CP@toff@wd=\CP@crop@wd
    \advance\CP@toff@wd     by -\CP@trim@wd
    \divide\CP@toff@wd  by \tw@
    \CP@toff@ht=\CP@crop@ht
    \advance\CP@toff@ht     by -\CP@trim@ht
    \divide\CP@toff@ht  by \tw@
    \advance\CP@trim@wd     by \CP@toff@wd
    \advance\CP@trim@ht     by \CP@toff@ht
    \CP@boff@wd=\CP@toff@wd
    \advance\CP@boff@wd     by -\CP@bled@wd
    \CP@boff@ht=\CP@toff@ht
    \advance\CP@boff@ht     by -\CP@bled@ht
    \advance\CP@bled@wd     by \CP@trim@wd
    \advance\CP@bled@ht     by \CP@trim@ht
}

\def\pdf@page@parameters{%
    \ifx\IsCropSet\str@yes
        [{ThisPage} << /CropBox [0 0 \strip@pt\CP@crop@wd\space \strip@pt\CP@crop@ht] >> /PUT pdfmark
        [{ThisPage} << /MediaBox[0 0 \strip@pt\CP@crop@wd\space \strip@pt\CP@crop@ht] >> /PUT pdfmark
    \fi
    \ifx\IsTrimSet\str@yes
        [{ThisPage} << /TrimBox [\strip@pt\CP@toff@wd\space \strip@pt\CP@toff@ht\space \strip@pt\CP@trim@wd\space \strip@pt\CP@trim@ht] >> /PUT pdfmark
    \fi
    \ifx\IsBleedSet\str@yes
        [{ThisPage} << /BleedBox[\strip@pt\CP@boff@wd\space \strip@pt\CP@boff@ht\space \strip@pt\CP@bled@wd\space \strip@pt\CP@bled@ht] >> /PUT pdfmark
    \fi
}

\def\shipout@PageObjects{%
    \special{ps: \pdf@page@parameters}%
}

\AtBeginDocument{%
    \do@pagebox@calc
    \let\org@begindvi\@begindvi
    \def\@begindvi{%
        \shipout@PageObjects
        \org@begindvi
        \global\let\@begindvi\rest@dvi@pages
    }%
}
\let\rest@dvi@pages\shipout@PageObjects

\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@}
\setlength\parindent{12pt}
\setlength\smallskipamount{3\p@ \@plus 1\p@ \@minus 1\p@}
\setlength\medskipamount{6\p@ \@plus 2\p@}
\setlength\bigskipamount{12\p@ \@plus 4\p@ \@minus 4\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\clubpenalty 10000
\widowpenalty 10000
\displaywidowpenalty 100
\predisplaypenalty   10000
\postdisplaypenalty  2500
\interlinepenalty 0
\brokenpenalty 10000
\lefthyphenmin=3
\righthyphenmin=3

% Contemporary large page layout
\setlength\headheight{16\p@}
\setlength\topmargin{2.24pc}
\addtolength\topmargin{-1in}
\setlength\topskip{10\p@}
\setlength\headsep{13.8\p@}
\setlength\footskip{16\p@}
\setlength\maxdepth{.5\topskip}
\setlength\textwidth{488.5pt}
\setlength\textheight{63\baselineskip}
\setlength\marginparsep{3\p@}
\setlength\marginparpush{3\p@}
\setlength\marginparwidth{35\p@}
\setlength\oddsidemargin{4.55pc}
\addtolength\oddsidemargin{-1in}
\setlength\@tempdima{\paperwidth}
\addtolength\@tempdima{-\textwidth}
\addtolength\@tempdima{-4.475pc}
\setlength\evensidemargin{\@tempdima}
\addtolength\evensidemargin{-1in}
\setlength\columnsep{18pt}
\setlength\columnseprule{0\p@}

\addtolength\textheight{\topskip}
\setlength\footnotesep{9\p@}
\setlength{\skip\footins}{12\p@ \@plus 6\p@ \@minus 1\p@}
\setcounter{totalnumber}{10}
\setcounter{topnumber}{5}
\setcounter{bottomnumber}{5}
\renewcommand\topfraction{.9}
\renewcommand\bottomfraction{.9}
\renewcommand\textfraction{.06}
\renewcommand\floatpagefraction{.94}
\renewcommand\dbltopfraction{.9}
\renewcommand\dblfloatpagefraction{.9}
\setlength\floatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\textfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}
\setlength\intextsep   {18\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dblfloatsep    {12\p@ \@plus 2\p@ \@minus 2\p@}
\setlength\dbltextfloatsep{20\p@ \@plus 2\p@ \@minus 4\p@}

\setlength\@fptop{0\p@}
\setlength\@fpsep{12\p@ \@plus 1fil}
\setlength\@fpbot{0\p@}
\setlength\@dblfptop{0\p@}
\setlength\@dblfpsep{12\p@ \@plus 1fil}
\setlength\@dblfpbot{0\p@}

\DeclareMathSizes{5}     {5}   {5}  {5}
\DeclareMathSizes{6}     {6}   {5}  {5}
\DeclareMathSizes{7}     {7}   {5}  {5}
\DeclareMathSizes{8}     {8}   {6}  {5}
\DeclareMathSizes{9}     {9}   {6.5}  {5}
\DeclareMathSizes{10}   {10}   {7.5}  {5}
\DeclareMathSizes{12}   {12}   {9}  {7}

\definecolor{jnlclr}{cmyk}{.80,.29,.05,0}
\definecolor{jnlruleclr}{cmyk}{.45,.06,.05,0}

% Contemporary page style
\def\ps@headings{%
  \let\@oddfoot\@empty
  \let\@evenfoot\@empty
  \def\@evenhead{\vbox{\hbox to \textwidth{\fontsize{7.8bp}{10}\selectfont\color{jnlclr}
        {\hspace*{1pc}{\selectfont\thepage}}\hspace*{8pt}\raisebox{-1.5pt}{\rule{.3pt}{8pt}}\hspace*{8pt}{\strut\leftmark}\hfill}\vspace{5\p@}{\color{black!15}\rule{\textwidth}{2\p@}}}}
  \def\@oddhead{\vbox{\hbox to \textwidth{\hfill\fontsize{7.8bp}{10}\selectfont\color{jnlclr}
        {{\strut\rightmark}}\hspace*{8pt}\raisebox{-1.5pt}{\rule{.3pt}{8pt}}\hspace*{8pt}{{\thepage}}\hspace*{1pc}}%
        \vspace{5\p@}{\color{black!15}\rule{\textwidth}{2\p@}}}}
  \def\titlemark##1{\markboth{##1}{##1}}
  \def\authormark##1{\gdef\leftmark{##1}}
}

\def\ps@opening{%
  \def\@oddfoot{{%
      \hbox to \textwidth{\parbox{\textwidth}{\hspace*{1pt}\vspace*{-0pt}\newline\rule{\textwidth}{.5pt}\vspace*{2pt}\newline%
      \fontsize{6.5}{8}\fontshape{n}\selectfont\raggedright \copyright\space The Author \@copyrightyear. Published by Oxford University Press. All rights reserved. For permissions, please e-mail: journals.permissions@oup.com%
        \hfill\vspace*{1pt}\newline\hspace*{1pt}\hfill\color{jnlclr}\textbf{\thepage}}}%
      }}
  \def\@evenfoot{{%
      \hbox to \textwidth{\parbox{\textwidth}{\hspace*{1pt}\vspace*{-0pt}\newline\rule{\textwidth}{.5pt}\vspace*{2pt}\newline%
      \fontsize{6.5}{8}\fontshape{n}\selectfont\raggedright \copyright\space The Author \@copyrightyear. Published by Oxford University Press. All rights reserved. For permissions, please e-mail: journals.permissions@oup.com%
        \hfill\vspace*{1pt}\newline\hspace*{1pt}\hfill\color{jnlclr}\textbf{\thepage}}}%
      }}
  \let\@evenhead\relax
  \let\@oddhead\relax
}

% Page range
\newif\iflastpagegiven   \lastpagegivenfalse
\newcommand\firstpage[1]{%
  \gdef\@firstpage{#1}%
  \ifnum\@firstpage>\c@page
    \setcounter{page}{#1}%
    \ClassWarning{BIO}{Increasing pagenumber to \@firstpage}%
  \else \ifnum\@firstpage<\c@page
    \ClassWarning{BIO}{Firstpage lower than pagenumber}\fi\fi
    \xdef\@firstpage{\the\c@page}%
    }
\def\@firstpage{1}
\def\pagenumbering#1{%
    \global\c@page \@ne
    \gdef\thepage{\csname @#1\endcsname \c@page}%
    \gdef\thefirstpage{\csname @#1\endcsname \@firstpage}%
    \gdef\thelastpage{\csname @#1\endcsname \@lastpage}%
    }

\newcommand\lastpage[1]{\xdef\@lastpage{#1}\global\lastpagegiventrue}
\def\@lastpage{0}
\def\setlastpage{\iflastpagegiven\else
    \edef\@tempa{@lastpage@}%
    \expandafter\ifx \csname \@tempa \endcsname \relax
        \gdef\@lastpage{0}%
    \else
        \xdef\@lastpage{\@nameuse{@lastpage@}}%
    \fi
    \fi}
\def\writelastpage{%
    \iflastpagegiven \else
    \immediate\write\@auxout%
    {\string\global\string\@namedef{@lastpage@}{\the\c@page}}%
    \fi}
\def\thepagerange{%
  \ifnum\@lastpage =0 {\ \bf ???} \else
  \ifnum\@lastpage = \@firstpage \ \thefirstpage\else
  \thefirstpage--\thelastpage \fi\fi}

\AtBeginDocument{\setlastpage\pagenumbering{arabic}}
\AtEndDocument{\writelastpage\if@final\clearemptydoublepage\else\clearpage\fi}

\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
\newcounter{subparagraph}[paragraph]
\newcounter{figure}
\newcounter{table}

\newenvironment{tablenotes}{\fontsize{6bp}{9}\selectfont\list{}{\setlength{\labelsep}{0pt}%
\setlength{\labelwidth}{0pt}%
\setlength{\leftmargin}{0pt}%
\setlength{\rightmargin}{0pt}%
\setlength{\topsep}{-6pt}%
\setlength{\itemsep}{2pt}%
\setlength{\partopsep}{0pt}%
\setlength{\listparindent}{0em}%
\setlength{\parsep}{0pt}}%
\item\relax%
}{\endlist\addvspace{0pt}}

\newcommand\thepage{\arabic{page}}
\renewcommand\thesection{\arabic{section}}
\renewcommand\thesubsection{{\thesection.\arabic{subsection}}}
\renewcommand\thesubsubsection{{\thesubsection.\arabic{subsubsection}}}
\renewcommand\theparagraph{\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph{\theparagraph.\arabic{subparagraph}}
\renewcommand\theequation{\arabic{equation}}

\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\newcommand\keywordsname{Key words:}
\newcommand\refname{References}
\newcommand\bibname{References}
\newcommand\indexname{Index}
\newcommand\figurename{Fig.}
\newcommand\tablename{Table}

\newcommand{\clearemptydoublepage}{\newpage{\pagestyle{empty}\cleardoublepage}}

\newif\if@mainmatter \@mainmattertrue
\newcommand\frontmatter{\clearpage\@mainmatterfalse\pagenumbering{roman}}
\newcommand\mainmatter{\clearpage\@mainmattertrue\pagenumbering{arabic}}
\newcommand\backmatter{\clearpage\@mainmatterfalse}

% Title
\newlength{\dropfromtop}
\setlength{\dropfromtop}{\z@}

\newif\if@appnotes
\newcommand{\application}{\global\@appnotestrue}

\long\def\title{\@ifnextchar[{\short@title}{\@@title}}
\def\short@title[#1]{\titlemark{#1}\@@@title}
\def\@@title#1{\authormark{#1}\@@@title{#1}}
\long\def\@@@title#1{\gdef\@title{#1}}
\def\@subtitle{}
\long\def\subtitle#1{\gdef\@subtitle{#1}}

\newcounter{myauthcount}
\setcounter{myauthcount}{0}
\regtotcounter{myauthcount}

\def\authorandsep{\ifnum\arabic{myauthcount@totc}=\arabic{myauthcount}\ifnum\arabic{myauthcount@totc}=1\else and \fi\else\fi}
\def\authorcommasep{\ifnum\arabic{myauthcount@totc}=\arabic{myauthcount}\else,\fi}

\def\author{\@@author}
\newcommand{\@@author}[2][]{%
  \g@addto@macro\@author{%
    \refstepcounter{myauthcount}%
      \hspace*{0.001pt}\hbox{\authorandsep#2\authorcommasep\ifx#1\@empty\else\textsuperscript{#1}\fi}\space}
}

\def\@vol{}
\def\vol#1{\global\def\@vol{#1}}
\def\@issue{}
\def\issue#1{\global\def\@issue{#1}}

\newcommand\defcase[1]{\@namedef{mycase@\the\numexpr#1\relax}}
\newcommand\myswitch[1]{\@nameuse{mycase@\the\numexpr#1\relax}}

\defcase{0}{Month}
\defcase{1}{January}
\defcase{2}{February}
\defcase{3}{March}
\defcase{4}{April}
\defcase{5}{May}
\defcase{6}{June}
\defcase{7}{July}
\defcase{8}{August}
\defcase{9}{September}
\defcase{10}{October}
\defcase{11}{November}
\defcase{12}{December}

\def\@history{}
\def\history#1{\global\def\@history{#1}}
\history{}

\newcounter{myhistorycount}
\setcounter{myhistorycount}{0}
\regtotcounter{myhistorycount}

\def\historycommasep{\ifnum\arabic{myhistorycount@totc}=\arabic{myhistorycount}\else;\ \fi}

\def\received#1#2#3{\g@addto@macro\@history{\refstepcounter{myhistorycount}FOR PUBLISHER ONLY Received on #1\ \myswitch{#2}\ #3\historycommasep}}
\def\revised#1#2#3{\g@addto@macro\@history{\refstepcounter{myhistorycount}revised on #1\ \myswitch{#2}\ #3\historycommasep}}
\def\accepted#1#2#3{\g@addto@macro\@history{\refstepcounter{myhistorycount}accepted on #1\ \myswitch{#2}\ #3\historycommasep}}

\def\abstract#1{\global\def\@abstract{#1}}

\def\keywords#1{\g@addto@macro\@abstract{%
    \vspace*{8pt}\par%
{\fontsize{8bp}{10}\textbf{{Key words}:}\ #1}%
}}

\def\@editor{}
\def\editor#1{\global\def\@editor{#1}}
\def\pubyear#1{\global\def\@pubyear{#1}}
\def\copyrightyear#1{\global\def\@copyrightyear{#1}}
\def\journaltitle#1{\global\def\@journaltitle{#1}}

\def\@boxedtext{}
\def\boxedtext{\@@boxedtext}
\newcommand{\@@boxedtext}[1]{\def\@boxedtext{\ifx#1\@empty\else\removelastskip\vskip-2pc\fbox{\parbox{.97\textwidth}{\textbf{Key Messages}\par#1\vspace{-8pt}}}\par\vskip2pc\fi}}

\def\orgdiv#1{#1}
\def\orgname#1{#1}
\def\orgaddress#1{#1}
\def\street#1{#1}
\def\postcode#1{#1}
\def\state#1{#1}
\def\country#1{#1}

\newcounter{myaddcount}
\setcounter{myaddcount}{0}
\regtotcounter{myaddcount}

\def\addressandsep{\ifnum\arabic{myaddcount@totc}=\arabic{myaddcount}\ifnum\arabic{myaddcount@totc}=1\else\ and \fi\else\fi}
\def\addresscommasep{\ifnum\arabic{myaddcount@totc}=\arabic{myaddcount}\else,\fi}

\def\@address{}
\def\address{\@@address}
\newcommand{\@@address}[2][]{%
  \g@addto@macro\@address{%
    \refstepcounter{myaddcount}%
    \ifnum\arabic{myaddcount}=1\vspace*{5pt}\par\fi%
      \addressandsep\ifx#1\@empty\else\textsuperscript{#1}\fi#2\addresscommasep}%
}
\def\@corresp{}
\def\corresp{\@@corresp}
\newcommand{\@@corresp}[2][]{%
  \g@addto@macro\@corresp{%
      \ifx#1\@empty\else\textsuperscript{#1}\fi#2}}

\def\@DOI{}
\def\DOI#1{\global\def\@DOI{#1}}

\definecolor{gray}{cmyk}{0, 0, 0, 0.15}
\definecolor{grayfifty}{cmyk}{0, 0, 0, 0.5}
\definecolor{graysixtyfive}{cmyk}{0, 0, 0, 0.65}
\newlength{\extraspace}
\setlength{\extraspace}{\z@}

\newcommand\maketitle{\par
  \begingroup
\def\@latex@warning@no@line##1{}%
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 3mm\noindent
        \@textsuperscript{\normalfont\@thefnmark}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage%
      \global\@topnum\z@
      \@maketitle
    \fi
    \thispagestyle{opening}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@address\@empty
  \global\let\@corresp\@empty
  \global\let\@history\@empty
  \global\let\@editor\@empty
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@subtitle\@empty
  \global\let\@title\@empty
  \global\let\@boxedtext\@empty
  \global\let\@pubyear\@empty
  \global\let\address\relax
  \global\let\boxedtext\relax
  \global\let\history\relax
  \global\let\editor\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\pubyear\relax
  \global\let\@copyrightline\@empty
  \global\let\and\relax
  \@afterindentfalse\@afterheading
}

\newlength{\aboveskipchk}
\setlength{\aboveskipchk}{\z@}

\def\@access{}
\def\access#1{\gdef\@access{#1}}

\def\appnotes#1{\gdef\@appnotes{#1}}
\appnotes{}

\newlength{\titlepagewd}
\setlength{\titlepagewd}{\textwidth}

\def\@maketitle{%
  \let\footnote\thanks
  \clearemptydoublepage
  \checkoddpage\ifcpoddpage\setlength{\aboveskipchk}{-28.6pt}\else\setlength{\aboveskipchk}{-28.6pt}\fi
  \vspace*{\aboveskipchk}%
  \vspace{\dropfromtop}%
  \vspace*{15.2pt}
  \if@twocolumn\else\vspace*{-10pt}\fi
  % Contemporary large header
  \hbox to \textwidth{{\color{black!20}\rule{45pt}{55pt}}\hfill
    \raisebox{-3.7pt}[0pt]{\hfill%
    \parbox[b]{416pt}{\raggedleft{{\fontsize{8bp}{10}\selectfont\itshape\hspace*{178pt} \@journaltitle,\rm\ \@copyrightyear, pp. \thepage--\thelastpage}}\\[-.5pt]
    \hspace*{181pt}{\color{black!15}\rule{184pt}{2pt}}\\[1.6pt]
    {\hfill\fontsize{6.5bp}{8}\selectfont\hspace*{178pt} \textbf{\color{jnlclr}doi: \@DOI}}\\[-.9pt]
    {\ifx\@access\@empty
    \else
    {\hfill\fontsize{6.5bp}{8}\selectfont\hspace*{178pt} \@access}\fi}
    \vskip-.9pt
    {\ifx\@appnotes\@empty
    \else
    {\hfill\fontsize{6.5bp}{8}\selectfont\hspace*{178pt} \@appnotes}\\[-.9pt]\fi}
      }}%
    }
  \vspace*{25pt}\par
  \hbox to \textwidth{%
  \parbox[t]{\titlepagewd}{%
    % Subtitle
    \ifx\@subtitle\@empty\else
    {\sffamilyfontcn\fontsize{14}{21}\selectfont\raggedright \@subtitle \par}%
    \vspace{7.5\p@}%
    \fi
    % Title
    {\fontsize{9.8bp}{10}\selectfont \MakeUppercase{\@appnotes}}\\[6.5pt]
    {\fontsize{18bp}{20}\bfseries\selectfont\leftskip0pt\rightskip0pt plus1fill {\color{jnlclr}\@title} \par}%
    \vspace{12\p@}%
    % Author
    {\boldmath\fontsize{13bp}{15}\selectfont\leftskip0pt\rightskip0pt plus1fill \@author \par}%
    \vspace{4\p@}%
    % Address
    {\fontsize{6.5bp}{9.5}\selectfont\leftskip0pt\rightskip0pt plus1fill \@address \par}%
    % Correspondence
    {\fontsize{6.5bp}{9.5}\selectfont\leftskip0pt\rightskip0pt plus1fill \@corresp \par}%
    % History
    {\fontsize{6.5bp}{9.5}\selectfont\leftskip0pt\rightskip0pt plus1fill \@history \par}%
    \vspace*{8.5pt}%
    % Abstract
    \vspace*{12.35pt}\par
    {\fontsize{11.5bp}{12}\fontseries{b}\selectfont\leftskip0pt\rightskip0pt plus1fill {\color{jnlclr} Abstract}\par}
    \vskip5.5pt
    \begingroup
    \begin{minipage}[t]{\textwidth}
    \parindent=0pt
    {\fontsize{9bp}{11}\selectfont
    \@abstract
    \vspace*{12pt}\par}
    \end{minipage}
    \endgroup
    }%
  }
  \vspace{13.5\p@}%
  \ifx\@boxedtext\@empty\else\vspace*{12pt}\par\fi%
  \@boxedtext
  \vspace{12\p@ plus 6\p@ minus 6\p@}%
  \vspace{\extraspace}
  \vspace*{-18pt}
}

% Abstract section
\newcommand{\absection}[1]{\par\noindent{\bfseries #1}\space\ignorespaces}

% Section macros
\setcounter{secnumdepth}{3}

\renewcommand{\@seccntformat}[1]{\if@unnumsec\else\csname the#1\endcsname.\space\fi}

% Contemporary section styles
\def\secsize{\fontsize{11bp}{13}\selectfont\fontseries{b}\leftskip0pt\rightskip0pt plus1fill\color{jnlclr}}
\def\subsecsize{\fontsize{10bp}{12}\selectfont\fontseries{b}\leftskip0pt\rightskip0pt plus1fill}
\def\subsubsecsize{\fontsize{9bp}{11}\selectfont\itshape\leftskip0pt\rightskip0pt plus1fill}

\def\section{%
  \@startsection{section}{1}{\z@}
  {-16\p@ plus -3\p@}{5\p@}
  {\reset@font\raggedright\secsize}}

\def\subsection{%
  \@startsection{subsection}{2}{\z@}
  {-11\p@ plus -2\p@}{2\p@}
  {\reset@font\raggedright\subsecsize}}

\def\subsubsection{%
  \@startsection{subsubsection}{3}{\z@}
  {-11\p@ plus -1\p@}{0.1em}
  {\reset@font\raggedright\subsubsecsize}}

\def\textcolon{\text{\rm :}}

\def\paragraph{%
  \@startsection{paragraph}{4}{\z@}
  {-6\p@}{0.1em}
  {\reset@font\normalsize}}

\def\@startsection#1#2#3#4#5#6{%
  \if@noskipsec \leavevmode \fi
  \par
  \@tempskipa #4\relax
  \@afterindenttrue
  \ifdim \@tempskipa <\z@
    \@tempskipa -\@tempskipa \@afterindentfalse
  \fi
  \if@nobreak
    \everypar{}%
  \else
    \addpenalty\@secpenalty\addvspace\@tempskipa
  \fi
  \@ifstar
    {\@ssect{#3}{#4}{#5}{#6}}%
    {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}}

\def\@sect#1#2#3#4#5#6[#7]#8{%
  \ifnum #2>\c@secnumdepth
    \let\@svsec\@empty
  \else
    \refstepcounter{#1}%
    \protected@edef\@svsec{\@seccntformat{#1}\relax}%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@
    \begingroup
      #6{%
        \@hangfrom{\hskip #3\relax\@svsec}%
          \interlinepenalty \@M #8\@@par}%
    \endgroup
    \csname #1mark\endcsname{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #7}%
  \else
    \def\@svsechd{%
      #6{\hskip #3\relax
      \@svsec #8}%
      \csname #1mark\endcsname{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #7}}%
  \fi
  \@xsect{#5}}

% Figures and tables
\setlength\arraycolsep{.5em}
\setlength\tabcolsep{.5em}
\setlength\arrayrulewidth{.5pt}
\setlength\doublerulesep{2.5pt}
\setlength\extrarowheight{\z@}
\renewcommand\arraystretch{1}

\setlength{\abovecaptionskip}{13pt}
\setlength{\belowcaptionskip}{2pt}

\long\def\@makecaption#1#2{\vspace{\abovecaptionskip}%
  \begingroup
  \scriptsize\sffamily
  \text{\sfb #1.}\space{#2}\par
  \endgroup}

\long\def\@tablecaption#1#2{%
  \begingroup%
  \fontsize{7pt}{9pt}\selectfont%
  \textbf{#1.}\space{#2\strut\par}%
  \endgroup\vspace{\belowcaptionskip}}

\long\def\@figurecaption#1#2{%
  \begingroup%
  \fontsize{6bp}{9pt}\selectfont%
  \hspace*{1pt}\vspace*{-3pt}\newline%
  \textbf{#1.}\space{#2\strut\par}%
  \endgroup\vspace{\belowcaptionskip}}

% Table rules
\def\toprule{\noalign{\ifnum0=`}\fi%
    \hrule \@height 0.5pt
    \hrule \@height 4pt \@width 0pt \futurelet \@tempa\@xhline}
\def\midrule{\noalign{\ifnum0=`}\fi%
    \hrule \@height 3pt \@width 0pt%
    \hrule \@height 0.5pt%
    \hrule \@height 4pt \@width 0pt \futurelet \@tempa\@xhline}
\def\botrule{\noalign{\ifnum0=`}\fi%
    \hrule \@height 3.75pt \@width 0pt
    \hrule \@height 0.5pt \futurelet \@tempa\@xhline}
\def\hrulefill{\leavevmode\leaders\hrule height .5pt\hfill\kern\z@}

\def\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\def\figure{\let\@makecaption\@figurecaption\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\let\@makecaption\@figurecaption\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}
\def\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{Table~\thetable}
\def\table{\let\@makecaption\@tablecaption\let\source\tablesource\@float{table}}
\def\endtable{\end@float}
\@namedef{table*}{\let\@makecaption\@tablecaption\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}

\newif\if@rotate \@rotatefalse
\newif\if@rotatecenter \@rotatecenterfalse
\def\rotatecenter{\global\@rotatecentertrue}
\def\rotateendcenter{\global\@rotatecenterfalse}
\def\rotate{\global\@rotatetrue}
\def\endrotate{\global\@rotatefalse}
\newdimen\rotdimen
\def\rotstart#1{\special{ps: gsave currentpoint currentpoint translate
    #1 neg exch neg exch translate}}
\def\rotfinish{\special{ps: currentpoint grestore moveto}}
\def\rotl#1{\rotdimen=\ht#1\advance\rotdimen by \dp#1
    \hbox to \rotdimen{\vbox to\wd#1{\vskip \wd#1
    \rotstart{270 rotate}\box #1\vss}\hss}\rotfinish}
\def\rotr#1{\rotdimen=\ht #1\advance\rotdimen by \dp#1
    \hbox to \rotdimen{\vbox to \wd#1{\vskip \wd#1
    \rotstart{90 rotate}\box #1\vss}\hss}\rotfinish}

\newdimen\tempdime
\newbox\temptbox

\begingroup
\catcode`\Q=3
\long\gdef\@ifmtarg#1{\@xifmtarg#1QQ\@secondoftwo\@firstoftwo\@nil}
\long\gdef\@xifmtarg#1#2Q#3#4#5\@nil{#4}
\long\gdef\@ifnotmtarg#1{\@xifmtarg#1QQ\@firstofone\@gobble\@nil}
\endgroup

\def\tablesize{\@setfontsize\tablesize{8\p@}{10\p@}}

\newenvironment{processtable}[3]{\setbox\temptbox=\hbox{{\tablesize #2}}%
\tempdime\wd\temptbox\@processtable{#1}{#2}{#3}{\tempdime}}
{\relax}

\newcommand{\@processtable}[4]{%
\if@rotate
\setbox4=\vbox to \hsize{\vss\hbox to \textheight{%
\begin{minipage}{#4}%
\@ifmtarg{#1}{}{\caption{#1}}{\tablesize #2}%
\vskip7\p@\noindent
\parbox{#4}{\fontsize{7}{9}\selectfont #3\par}%
\end{minipage}}\vss}%
\rotr{4}
\else
\hbox to \hsize{\hss\begin{minipage}[t]{#4}%
\vskip2.9pt
\@ifmtarg{#1}{}{\caption{#1}}{\tablesize #2}%
\vskip6\p@\parindent=12pt
\parbox{#4}{\fontsize{7}{9}\selectfont #3\par}%
\end{minipage}\hss}\fi}

\newcolumntype{P}[1]{>{\raggedright\let\\\@arraycr\hangindent1em}p{#1}}

% Lists
\def\labelenumi{{\rm\arabic{enumi}.}}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{{\rm\alph{enumii}.}}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{{\rm(\roman{enumiii})}}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{{\rm(\arabic{enumiv})}}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}
\def\labelitemi{{\small$\bullet$}}
\def\labelitemii{{--}}
\def\labelitemiii{{\small$\bullet$}}
\def\labelitemiv{{\small$\bullet$}}

\def\@listI{\leftmargin\leftmargini \topsep\medskipamount}
\let\@listi\@listI
\@listi
\def\@listii{\topsep\z@\leftmargin\leftmarginii}
\def\@listiii{\leftmargin\leftmarginiii \topsep\z@}
\def\@listiv{\leftmargin\leftmarginiv \topsep\z@}
\def\@listv{\leftmargin\leftmarginv \topsep\z@}
\def\@listvi{\leftmargin\leftmarginvi \topsep\z@}

\setlength{\leftmargini}{3mm}
\setlength{\leftmarginii}{1mm}
\setlength{\leftmarginiii}{2mm}
\setlength{\leftmarginiv}{\z@}

\def\enumargs{%
   \partopsep     \z@
   \itemsep       \z@
   \parsep        \z@
   \labelsep      1em
   \listparindent \parindent
   \itemindent    \z@
   \topsep 7\p@
}

\def\enumerate{\@ifnextchar[{\@numerate}{\@numerate[0]}}

\def\@numerate[#1]{%
     \ifnum \@enumdepth >3 \@toodeep\else
     \advance\@enumdepth \@ne
     \edef\@enumctr{enum\romannumeral\the\@enumdepth}
     \list{\csname label\@enumctr\endcsname}{%
       \enumargs
       \setlength{\leftmargin}{\csname leftmargin\romannumeral\the\@enumdepth\endcsname}
       \usecounter{\@enumctr}
       \settowidth\labelwidth{#1}
       \addtolength{\leftmargin}{\labelwidth}
       \addtolength{\leftmargin}{2pt}
       \def\makelabel##1{\hss \llap{##1}}}%
     \fi
   }
\let\endenumerate\endlist

\def\itemargs{%
   \partopsep     \z@
   \itemsep       0\p@
   \parsep        \z@
   \labelsep      1em
   \rightmargin   \z@
   \listparindent \parindent
   \itemindent    \z@
   \topsep7\p@
   \setlength{\leftmarginii}{3mm}
}

\def\itemize{\@ifnextchar[{\@itemize}{\@itemize[$\bullet$]}}

\def\@itemize[#1]{%
     \ifnum \@itemdepth >3 \@toodeep\else
     \advance\@itemdepth \@ne
     \edef\@itemctr{item\romannumeral\the\@itemdepth}
     \list{\csname label\@itemctr\endcsname}{%
       \itemargs
       \setlength{\leftmargin}{\csname leftmargin\romannumeral\the\@itemdepth\endcsname}
       \settowidth\labelwidth{#1}
       \addtolength{\leftmargin}{\labelwidth}
       \def\makelabel##1{\hss \llap{##1}}}%
     \fi
   }
\let\enditemize\endlist

\newenvironment{unlist}{%
  \begin{list}{}%
    {\setlength{\labelwidth}{\z@}%
     \setlength{\labelsep}{\z@}%
     \setlength{\topsep}{\medskipamount}%
     \setlength{\itemsep}{3\p@}%
     \setlength{\leftmargin}{2em}%
     \setlength{\itemindent}{-2em}}}
{\end{list}}

% Quotes
\def\quotation{\par\begin{list}{}{
      \setlength{\topsep}{\medskipamount}
      \setlength{\leftmargin}{2em}%
      \setlength{\rightmargin}{\z@}%
      \setlength\labelwidth{0pt}%
      \setlength\labelsep{0pt}%
      \listparindent\parindent}%
  \item[]}
\def\endquotation{\end{list}}
\let\quote\quotation
\let\endquote\endquotation

\skip\@mpfootins = \skip\footins
\fboxsep=6\p@
\fboxrule=1\p@

% Table of contents
\newcommand\@pnumwidth{4em}
\newcommand\@tocrmarg{2.55em plus 1fil}
\newcommand\@dotsep{1000}
\setcounter{tocdepth}{4}

\def\numberline#1{\hbox to \@tempdima{{#1}}}

\def\@authortocline#1#2#3#4#5{%
  \vskip 1.5\p@
  \ifnum #1>\c@tocdepth \else
    {\leftskip #2\relax \rightskip \@tocrmarg \parfillskip -\rightskip
     \parindent #2\relax\@afterindenttrue
     \interlinepenalty\@M
     \leavevmode
     \@tempdima #3\relax
     \advance\leftskip \@tempdima \null\nobreak\hskip -\leftskip
     {\itshape #4}\nobreak
     \leaders\hbox{$\m@th
        \mkern \@dotsep mu\hbox{.}\mkern \@dotsep
        mu$}\hfill
     \nobreak
     \hb@xt@\@pnumwidth{\hfil}%
     \par}%
  \fi}

\newcommand*\l@author{\@authortocline{2}{0pt}{30pt}}
\newcommand*\l@section{\@dottedtocline{3}{11pt}{20pt}}
\newcommand*\l@subsection{\@dottedtocline{4}{31pt}{29pt}}
\newcommand*\l@subsubsection[2]{}

% Footnotes
\def\footnoterule{\noindent\rule{\columnwidth}{0.5pt}}
\def\@makefnmark{\@textsuperscript{\normalfont\@thefnmark}}
\newcommand\@makefntext[1]{\noindent{\@makefnmark\enskip}#1}

% References
\providecommand{\newblock}{}
\newenvironment{thebibliography}{%
  \section{\bibname}%
  \begingroup
  \small
  \begin{list}{}{%
      \setlength{\topsep}{\z@}%
      \setlength{\labelsep}{\z@}%
      \settowidth{\labelwidth}{\z@}%
      \setlength{\leftmargin}{4mm}%
      \setlength{\itemindent}{-4mm}}\small}
{\end{list}\endgroup}

\RequirePackage{natbib}

\def\NAT@force@numbers{\if@numbib\global\NAT@numberstrue\else\global\NAT@numbersfalse\fi}

\renewcommand\@biblabel[1]{#1.\hfill\hspace*{3pt}}

\renewcommand\NAT@bibsetnum[1]{\settowidth\labelwidth{\@biblabel{#1}}%
   \setlength{\leftmargin}{\labelwidth}\addtolength{\leftmargin}{\labelsep}%
   \setlength{\itemsep}{\bibsep}\setlength{\parsep}{\z@}%
   \ifNAT@openbib
     \addtolength{\leftmargin}{\bibindent}%
     \setlength{\itemindent}{-\bibindent}%
     \setlength{\listparindent}{\itemindent}%
     \setlength{\parsep}{0pt}%
   \fi
}

% Math parameters
\setlength{\jot}{5\p@}
\mathchardef\@m=1500

\def\frenchspacing{\sfcode`\.\@m \sfcode`\?\@m \sfcode`\!\@m
  \sfcode`\:\@m \sfcode`\;\@m \sfcode`\,\@m}

% Theorems
\def\th@plain{%
\thm@headfont{\quad\scshape}
\thm@notefont{\upshape\mdseries}
\thm@headpunct{.}
\thm@headsep 5\p@ plus\p@ minus\p@\relax
\itshape
}

\vbadness=9999
\tolerance=9999
\doublehyphendemerits=10000
\doublehyphendemerits 640000
\finalhyphendemerits  1000000

\flushbottom
\frenchspacing
\ps@headings
\twocolumn

\newcommand{\medline}[1]{\unskip\unskip\ignorespaces}

% Methods environment
\newenvironment{methods}{%
  \begingroup
\def\section{%
  \@startsection{section}{1}{\z@}
  {-24\p@ plus -3\p@}{4\p@}
  {\reset@font\raggedright\sffamilyfontbold\fontsize{10}{12}\selectfont}}
 \def\subsection{%
  \@startsection{subsection}{2}{\z@}
  {-11\p@ plus -2\p@}{4\p@}
  {\reset@font\raggedright\sffamilyfont\fontsize{9}{12}\selectfont}}
\def\subsubsection{%
  \@startsection{subsubsection}{3}{\z@}
  {-11\p@ plus -1\p@}{0.001em}
  {\reset@font\normalfont\mathversion{bold}\normalsize\bfseries}}
\normalsize
  \par}
{\par\endgroup\bigskip\@afterheading\@afterindentfalse}

\language=2
\hyphenation{Figure Table Figures Tables}

% Biography
\RequirePackage{wrapfig}

\newcount\wraplines
\wraplines=8

\newbox\@authorfigbox
\newskip\@authorfigboxdim
\newskip\biofigadjskip
\biofigadjskip=0pt

\def\authbiotextfont{\reset@font\fontsize{8bp}{9.5bp}\selectfont}

\newenvironment{biography}[2]{\par\addvspace{11.5pt plus3.375pt minus1.6875pt}%
\def\author##1{{\bfseries##1}}%
\if!#1!\def\@authorfig{}\else\def\@authorfig{{#1}}\fi%
\setbox\@authorfigbox=\hbox{#1}%
\@authorfigboxdim=\wd\@authorfigbox%
\advance\@authorfigboxdim by -10pt
\wraplines=9\fboxrule=1pt\fboxsep=6pt%
\noindent{%
\ifx\@authorfig\@empty\else\unskip%
\begin{wrapfigure}[\wraplines]{l}[0pt]{\@authorfigboxdim}%
\vskip-19pt\vskip\biofigadjskip%
\@authorfig%
\end{wrapfigure}%
\fi%
{\authbiotextfont#2\par}%
\par%
}}{\par\addvspace{10.5pt plus3.375pt minus1.6875pt}}

% Theorem environments
\RequirePackage{amsthm}

\@ifpackageloaded{amsthm}{%
\let\proof\relax
\let\endproof\relax

\def\@begintheorem#1#2[#3]{%
  \deferred@thm@head{\the\thm@headfont \thm@indent
    \@ifempty{#1}{\let\thmname\@gobble}{\let\thmname\@iden}%
    \@ifempty{#2}{\let\thmnumber\@gobble}{\let\thmnumber\@iden}%
    \@ifempty{#3}{\let\thmnote\@gobble}{\let\thmnote\@iden}%
    \thm@swap\swappedhead\thmhead{#1}{#2}{#3}%
    \the\thm@headpunct
    \thmheadnl
    \hskip\thm@headsep
  }%
  \ignorespaces
}

\def\@endtheorem{\endtrivlist\@endpefalse}

\AtBeginDocument{%
\DeclareSymbolFont{AMSa}{U}{msa}{m}{n}%
\DeclareMathSymbol{\opensquare}{\mathord}{AMSa}{"03}%
\def\qedsymbol{\ensuremath{\opensquare}}%
\newenvironment{proof}[1][\proofname]{\par\removelastskip%
  \pushQED{\ensuremath{\qed}}%
  \normalfont \topsep7.5\p@\@plus7.5\p@\relax%
  \trivlist%
  \item[\hskip\labelsep%
        \itshape%
    #1\ \@addpunct{}]\ignorespaces%
}{%
  \popQED\endtrivlist\@endpefalse%
}}

\def\thm@space@setup{%
\thm@preskip=12pt%
\thm@postskip=12pt}

\newtheoremstyle{thmstyleone}%
{12pt plus2pt minus1pt}{12pt plus2pt minus1pt}%
{\normalfont\itshape}{0pt}{\bfseries}{}{.5em}%
{\thmname{#1}\thmnumber{\@ifnotempty{#1}{ }\@upn{#2}}\thmnote{ {\the\thm@notefont(#3)}}}

\newtheoremstyle{thmstyletwo}%
{12pt plus2pt minus1pt}{12pt plus2pt minus1pt}%
{\itshape}{0pt}{\normalfont}{}{.5em}%
{\thmname{#1}\thmnumber{\@ifnotempty{#1}{ }{#2}}\thmnote{ {\the\thm@notefont(#3)}}}

\newtheoremstyle{thmstylethree}%
{12pt plus2pt minus1pt}{12pt plus2pt minus1pt}%
{\normalfont}{0pt}{\bfseries}{}{.5em}%
{\thmname{#1}\thmnumber{\@ifnotempty{#1}{ }\@upn{#2}}\thmnote{ {\the\thm@notefont(#3)}}}

\newtheoremstyle{thmstylefour}%
{12pt plus2pt minus1pt}{12pt plus2pt minus1pt}%
{\normalfont}{0pt}{\itshape}{}{.5em}%
{\global\proofthmtrue\thmname{#1} \thmnote{#3}}
}{}

\def\sbond{\ensuremath{\raise.25ex\hbox{${-}\!\!\!\!{-}$}}\kern -.9pt}
\def\dbond{\ensuremath{\raise.25ex\hbox{=$\!$=}}}
\def\tbond{\ensuremath{\raise.20ex\hbox{${\equiv}\!\!\!{\equiv}$}}}

\newcommand{\query}[2][0pt]{}

\renewcommand{\dag}{{\mathversion{normal}$^{\dagger}$}}

% Additional Packages
\RequirePackage[countmax]{subfloat}
\usepackage{anyfontsize}
\RequirePackage{multirow}
\RequirePackage{footnote}
\RequirePackage{url}
\RequirePackage{amsmath}
\RequirePackage{mathrsfs}
\RequirePackage{algorithm}
\RequirePackage{algorithmicx}
\RequirePackage{algpseudocode}
\RequirePackage{listings}

% Appendix macros
\newif\if@chapter@pp\@chapter@ppfalse
\newif\if@knownclass@pp\@knownclass@ppfalse
\@ifundefined{chapter}{%
  \@ifundefined{section}{}{\@knownclass@pptrue}}{%
  \@chapter@pptrue\@knownclass@pptrue}
\providecommand{\phantomsection}{}
\newcounter{@pps}
  \renewcommand{\the@pps}{\alph{@pps}}
\newif\if@pphyper
  \@pphyperfalse
\AtBeginDocument{\@ifpackageloaded{hyperref}{\@pphypertrue}{}}

\newif\if@dotoc@pp\@dotoc@ppfalse
\newif\if@dotitle@pp\@dotitle@ppfalse
\newif\if@dotitletoc@pp\@dotitletoc@ppfalse
\newif\if@dohead@pp\@dohead@ppfalse
\newif\if@dopage@pp\@dopage@ppfalse

\newcommand{\@ppendinput}{}
\if@knownclass@pp\else
  \renewcommand{\@ppendinput}{\endinput}
\fi
\@ppendinput

\newcommand{\appendixtocon}{\@dotoc@pptrue}
\newcommand{\appendixtocoff}{\@dotoc@ppfalse}
\newcommand{\appendixpageon}{\@dopage@pptrue}
\newcommand{\appendixpageoff}{\@dopage@ppfalse}
\newcommand{\appendixtitleon}{\@dotitle@pptrue}
\newcommand{\appendixtitleoff}{\@dotitle@ppfalse}
\newcommand{\appendixtitletocon}{\@dotitletoc@pptrue}
\newcommand{\appendixtitletocoff}{\@dotitletoc@ppfalse}
\newcommand{\appendixheaderon}{\@dohead@pptrue}
\newcommand{\appendixheaderoff}{\@dohead@ppfalse}
\newcounter{@ppsavesec}
\newcounter{@ppsaveapp}
\setcounter{@ppsaveapp}{0}
\newcommand{\@ppsavesec}{%
  \if@chapter@pp \setcounter{@ppsavesec}{\value{chapter}} \else
                 \setcounter{@ppsavesec}{\value{section}} \fi}
\newcommand{\@pprestoresec}{%
  \if@chapter@pp \setcounter{chapter}{\value{@ppsavesec}} \else
                 \setcounter{section}{\value{@ppsavesec}} \fi}
\newcommand{\@ppsaveapp}{%
  \if@chapter@pp \setcounter{@ppsaveapp}{\value{chapter}} \else
                 \setcounter{@ppsaveapp}{\value{section}} \fi}
\newcommand{\restoreapp}{%
  \if@chapter@pp \setcounter{chapter}{\value{@ppsaveapp}} \else
                 \setcounter{section}{\value{@ppsaveapp}} \fi}
\providecommand{\appendixname}{Appendix}
\newcommand{\appendixtocname}{Appendices}
\newcommand{\appendixpagename}{Appendices}
\newcommand{\appendixpage}{\if@chapter@pp \@chap@pppage \else \@sec@pppage \fi}
\newcommand{\clear@ppage}{\if@openright\cleardoublepage\else\clearpage\fi}

\newcommand{\@chap@pppage}{%
  \clear@ppage
  \thispagestyle{plain}%
  \if@twocolumn\onecolumn\@tempswatrue\else\@tempswafalse\fi
  \null\vfil
  \markboth{}{}%
  {\centering
   \interlinepenalty \@M
   \normalfont
   \Huge \bfseries \appendixpagename\par}%
  \if@dotoc@pp
    \addappheadtotoc
  \fi
  \vfil\newpage
  \if@twoside
    \if@openright
      \null
      \thispagestyle{empty}%
      \newpage
    \fi
  \fi
  \if@tempswa
    \twocolumn
  \fi
}

\newcommand{\@sec@pppage}{%
  \par
  \addvspace{4ex}%
  \@afterindentfalse
  {\parindent \z@ \raggedright
   \interlinepenalty \@M
   \normalfont
   \huge \bfseries \appendixpagename%
   \markboth{}{}\par}%
  \if@dotoc@pp
    \addappheadtotoc
  \fi
  \nobreak
  \vskip 3ex
  \@afterheading
}

\newif\if@pptocpage
  \@pptocpagetrue
\newcommand{\noappendicestocpagenum}{\@pptocpagefalse}
\newcommand{\appendicestocpagenum}{\@pptocpagetrue}
\newcommand{\addappheadtotoc}{%
  \phantomsection
  \if@chapter@pp
    \if@pptocpage
      \addcontentsline{toc}{chapter}{\appendixtocname}%
    \else
      \if@pphyper
        \addtocontents{toc}{\protect\contentsline{chapter}{\appendixtocname}{}{\@currentHref}}%
      \else
        \addtocontents{toc}{\protect\contentsline{chapter}{\appendixtocname}{}}%
      \fi
    \fi
  \else
    \if@pptocpage
      \addcontentsline{toc}{section}{\appendixtocname}%
    \else
      \if@pphyper
        \addtocontents{toc}{\protect\contentsline{section}{\appendixtocname}{}{\@currentHref}}%
      \else
        \addtocontents{toc}{\protect\contentsline{section}{\appendixtocname}{}}%
      \fi
    \fi
  \fi
}

\providecommand{\theH@pps}{\alph{@pps}}

\newcommand{\@resets@pp}{\par
  \@ppsavesec
  \stepcounter{@pps}
  \setcounter{section}{0}%
  \if@chapter@pp
    \setcounter{chapter}{0}%
    \renewcommand\@chapapp{\appendixname}%
    \renewcommand\thechapter{\@Alph\c@chapter}%
  \else
    \setcounter{subsection}{0}%
    \renewcommand\thesection{\@Alph\c@section}%
  \fi
  \if@pphyper
    \if@chapter@pp
      \renewcommand{\theHchapter}{\theH@pps.\Alph{chapter}}%
    \else
      \renewcommand{\theHsection}{\theH@pps.\Alph{section}}%
    \fi
    \xdef\Hy@chapapp{\Hy@appendixstring}%
  \fi
  \restoreapp
}

\newenvironment{appendices}{%
  \@resets@pp
  \if@dotoc@pp
    \if@dopage@pp
      \if@chapter@pp
        \clear@ppage
      \fi
      \appendixpage
    \else
       \if@chapter@pp
         \clear@ppage
       \fi
      \addappheadtotoc
    \fi
  \else
    \if@dopage@pp
      \appendixpage
    \fi
  \fi
  \if@chapter@pp
    \if@dotitletoc@pp \@redotocentry@pp{chapter} \fi
  \else
    \if@dotitletoc@pp \@redotocentry@pp{section} \fi
    \if@dohead@pp
      \def\sectionmark##1{%
        \if@twoside
          \markboth{\@formatsecmark@pp{##1}}{}
        \else
          \markright{\@formatsecmark@pp{##1}}{}
        \fi}
    \fi
    \if@dotitle@pp
      \def\sectionname{\appendixname}
      \def\@seccntformat##1{\@ifundefined{##1name}{}{\csname ##1name\endcsname\ }%
        \csname the##1\endcsname\quad}
    \fi
  \fi}{%
  \@ppsaveapp\@pprestoresec}

\newcommand{\setthesection}{\thechapter.\Alph{section}}
\newcommand{\setthesubsection}{\thesection.\Alph{subsection}}

\newcommand{\@resets@ppsub}{\par
  \stepcounter{@pps}
  \if@chapter@pp
    \setcounter{section}{0}
    \renewcommand{\thesection}{\setthesection}
  \else
    \setcounter{subsection}{0}
    \renewcommand{\thesubsection}{\setthesubsection}
  \fi
  \if@pphyper
    \if@chapter@pp
      \renewcommand{\theHsection}{\theH@pps.\setthesection}%
    \else
      \renewcommand{\theHsubsection}{\theH@pps.\setthesubsection}%
    \fi
    \xdef\Hy@chapapp{\Hy@appendixstring}%
  \fi
}

\newenvironment{subappendices}{%
  \@resets@ppsub
  \if@chapter@pp
    \if@dotitletoc@pp \@redotocentry@pp{section} \fi
    \if@dotitle@pp
      \def\sectionname{\appendixname}
      \def\@seccntformat##1{\@ifundefined{##1name}{}{\csname ##1name\endcsname\ }%
        \csname the##1\endcsname\quad}
    \fi
  \else
    \if@dotitletoc@pp \@redotocentry@pp{subsection} \fi
    \if@dotitle@pp
      \def\subsectionname{\appendixname}
      \def\@seccntformat##1{\@ifundefined{##1name}{}{\csname ##1name\endcsname\ }%
        \csname the##1\endcsname\quad}
    \fi
  \fi}{}

\newcommand{\@formatsecmark@pp}[1]{%
  \MakeUppercase{\appendixname\space
    \ifnum \c@secnumdepth >\z@
      \thesection\quad
    \fi
    #1}}
\newcommand{\@redotocentry@pp}[1]{%
  \let\oldacl@pp=\addcontentsline
  \def\addcontentsline##1##2##3{%
    \def\@pptempa{##1}\def\@pptempb{toc}%
    \ifx\@pptempa\@pptempb
      \def\@pptempa{##2}\def\@pptempb{#1}%
      \ifx\@pptempa\@pptempb
\oldacl@pp{##1}{##2}{\appendixname\space ##3}%
      \else
        \oldacl@pp{##1}{##2}{##3}%
      \fi
    \else
      \oldacl@pp{##1}{##2}{##3}%
    \fi}
}

% Hyperref
\RequirePackage[colorlinks,linkcolor=blue,urlcolor=blue,citecolor=blue,anchorcolor=blue]{hyperref}

% ORCID logo
\RequirePackage{tikz}
\usetikzlibrary{svg.path}

\definecolor{orcidlogocol}{HTML}{A6CE39}
\tikzset{
  orcidlogo/.pic={
    \fill[orcidlogocol] svg{M256,128c0,70.7-57.3,128-128,128C57.3,256,0,198.7,0,128C0,57.3,57.3,0,128,0C198.7,0,256,57.3,256,128z};
    \fill[white] svg{M86.3,186.2H70.9V79.1h15.4v48.4V186.2z}
                 svg{M108.9,79.1h41.6c39.6,0,57,28.3,57,53.6c0,27.5-21.5,53.6-56.8,53.6h-41.8V79.1z M124.3,172.4h24.5c34.9,0,42.9-26.5,42.9-39.7c0-21.5-13.7-39.7-43.7-39.7h-23.7V172.4z}
                 svg{M88.7,56.8c0,5.5-4.5,10.1-10.1,10.1c-5.6,0-10.1-4.6-10.1-10.1c0-5.6,4.5-10.1,10.1-10.1C84.2,46.7,88.7,51.3,88.7,56.8z};
  }
}

\newcommand{\@OrigHeightRecip}{0.00390625}
\newlength{\@curXheight}

\DeclareRobustCommand\ORCID[1]{%
\texorpdfstring{%
\setlength{\@curXheight}{\fontcharht\font`X}%
\href{https://orcid.org/#1}{\XeTeXLinkBox{\mbox{%
\begin{tikzpicture}[yscale=-\@OrigHeightRecip*\@curXheight,
xscale=\@OrigHeightRecip*\@curXheight,transform shape]
\pic{orcidlogo};
\end{tikzpicture}%
}}}}{}}

\def \@fpsadddefault {%
   \edef \@fps {\@fps\csname fps@\@captype \endcsname}%
}

\endinput
