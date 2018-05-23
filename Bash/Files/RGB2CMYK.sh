#!/bin/bash

## Convert RGB pdfs to CMYK
gs -dSAFER -dBATCH -dNOPAUSE -dNOCACHE \
	-sDEVICE=pdfwrite \
	-sColorConversionStrategy=CMYK \
	-sColorConversionStrategyForImages=CMYK \
	-dProcessColorModel=/DeviceCMYK \
	-dPDFSETTINGS=/prepress \
	-sOutputFile="${1/.pdf}"-CMYK.pdf "$1"
