#! /bin/bash

ROOTFILE="$1"
PDFFILE="$2"

MACRO="muonGlobalAlignmentMftDCAtest.C"

WD=$(dirname $0)

if [ -z "$ROOTFILE" ]; then
    root -l -b -q "$WD/$MACRO()"
elif [ -z "$PDFFILE" ]; then
    root -l -b -q "$WD/$MACRO(\"$ROOTFILE\")"
else
    root -l -b -q "$WD/$MACRO(\"$ROOTFILE\", \"$PDFFILE\")"
fi
