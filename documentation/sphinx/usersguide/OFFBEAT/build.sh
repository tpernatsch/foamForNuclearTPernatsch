#!/bin/bash
cd ${0%/*} || exit 1    # run from this directory

SRC=../Cases/Verification
TARGET=docs/test_cases/Verification

echo Copying test case results to the docs folder
rm -rf $TARGET
mkdir -p $TARGET
(cd $SRC && tar -cf - $(find . -type f | egrep "(.svg|.jpg|.png|.md)$")) | (cd $TARGET && tar -x)

echo Generating docs
python generate_docs.py

echo Building docs using mkdocs
mkdocs build
