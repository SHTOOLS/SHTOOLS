!#/bin/bash

echo 'Killing all Jekyll instances'
kill -9 $(ps aux | grep '[j]ekyll' | awk '{print $2}')
clear

echo "Building PDF-friendly HTML site for Mydoc ...";
bundle exec jekyll serve --detach --config _config.yml,pdfconfigs/config_mydoc_pdf.yml;
echo "done";

echo "Building the PDF ...";
prince --javascript --http-timeout=10 --input-list=_site/pdfconfigs/prince-list.txt -o pdf/shtools.pdf;

echo "Done. Look in the pdf directory to see if it printed successfully."
