#!/bin/sh

# obtain a recent ClinVar release, first grab the XML dump
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_2014-04.xml.gz

# obtain a smaller tab-separated ClinVar dump, not sure whether I will use the 
# tab-separated dump, or the XML file yet
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
