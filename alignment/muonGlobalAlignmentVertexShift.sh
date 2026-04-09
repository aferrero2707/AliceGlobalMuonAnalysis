#! /bin/bash

ROOTFILE="$1"
PDFFILE="$2"

WD=$(dirname $0)

if [ -z "$ROOTFILE" ]; then
    root -l -b -q "$WD/muonGlobalAlignmentVertexShift.C()"
elif [ -z "$PDFFILE" ]; then
    root -l -b -q "$WD/muonGlobalAlignmentVertexShift.C(\"$ROOTFILE\")"
else
    root -l -b -q "$WD/muonGlobalAlignmentVertexShift.C(\"$ROOTFILE\", \"$PDFFILE\")"
fi
