#! /bin/bash

DE=$1
ROOTFILE="$2"

WD=$(dirname $0)

if [ -z "$DE" ]; then
    echo "DE not specified, exiting"
    exit 1
elif [ -z "$ROOTFILE" ]; then
    root -l -b -q "$WD/muonGlobalAlignmentDEResiduals.C($DE)"
else
    root -l -b -q "$WD/muonGlobalAlignmentDEResiduals.C($DE, \"$ROOTFILE\")"
fi
