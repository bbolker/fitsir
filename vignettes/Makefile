%.pdf: %.Rnw
	echo "knitr::knit2pdf(\"$<\")" | R --slave

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.tex *.o *.so
