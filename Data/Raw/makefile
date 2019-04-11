export PATH := /Library/TeX/texbin:$(PATH)    # add LaTeX path
export PATH := /usr/local/bin:~/.local/bin:$(PATH) # add pandoc-citeproc-preamble

# Cluster targets
pdf:  dataReportTemplate.pdf
html:  dataReportTemplate.html

all: pdf html


# # Variables for pandoc calls
# PANDOC = pandoc
# LATEX = pdflatex
# BIB  = --bibliography=/Users/sdnjohnson/Work/Library/library.bib 

# CSL  = --csl=/Users/sdnjohnson/Work/write/templates/csl/cjfas.csl 
# # CSL repository https://github.com/citation-style-language/styles

# PTMP = --template=/Users/sdnjohnson/Work/write/templates/tex/basic.latex
# NOPRETMP = --template=/Users/sdnjohnson/Work/write/templates/tex/noPreamble.latex
# # BIBPRE = --filter pandoc-citeproc-preamble -M citeproc-preamble=/Users/sdnjohnson/Work/write/templates/tex/citeproc-preamble.latex 
# EQNO = --filter pandoc-eqnos

# ABV  = --citation-abbreviations=/Users/sdnjohnson/Work/write/templates/csl/csl-abbrev/jabbrev.json
# # STMP = --template ./pandoc/sup.latex
# # FTMP = --template ./pandoc/fig.latex
# FONT = --variable fontsize=12pt
# TIMS = --variable mainfont=Times
# WTMP = --reference-doc=/Users/sdnjohnson/Work/write/templates/docx/basic.docx

# Can create word doc with formatting want and save for template for this

## pdf
%.pdf: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", clean = FALSE, output_format = "bookdown::pdf_document2")'

# html
%.html: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", output_format = "bookdown::html_document2")'

# docx
%.docx: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", output_format = "word_document")'

# all formats (depends on yaml header)
%.all: %.Rmd makefile
	Rscript -e 'rmarkdown::render(input = "$<", output_format = "all")'


## appendix.tex
# appendix.tex: appendix.md makefile
# 	$(PANDOC) -s $(NOPRETMP) $(ABV) $(BIB) $(BIBPRE) $(CSL) $(FONT) $(TIMS) -o $@ $<

# TandF
# TandF.pdf: TandF.tex makefile
# 	$(LATEX) $<

# TandF.tex: TandF.md makefile
# 	$(PANDOC) -s -S $(PTMP) -o $@ $<	

# TandFnopre.tex: TandF.md makefile
# 	$(PANDOC) -s -S $(NOPRETMP) -o $@ $<	

## cleanup
clean:
	rm -f *.{pdf,doc,docx,tex,aux,log,out,toc,bbl}
