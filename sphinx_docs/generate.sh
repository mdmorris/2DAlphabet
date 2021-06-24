rm source/doctrees/*.rst
sphinx-apidoc -o source/doctrees/ -H="2DAlphabet API" -f -d 3 -e -P --templatedir source/_templates/ ../TwoDAlphabet
make html