#! /bin/bash

ROOTFILE="$1"
PDFFILE="$2"

WD=$(dirname $0)

if [ -z "$ROOTFILE" ]; then
    root -l -b -q "$WD/muonGlobalAlignmentMftMchResiduals.C()"
elif [ -z "$PDFFILE" ]; then
    root -l -b -q "$WD/muonGlobalAlignmentMftMchResiduals.C(\"$ROOTFILE\")"
else
    root -l -b -q "$WD/muonGlobalAlignmentMftMchResiduals.C(\"$ROOTFILE\", \"$PDFFILE\")"
fi
