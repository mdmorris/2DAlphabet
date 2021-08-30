rm source/doctrees/*.rst
sphinx-apidoc -o source/doctrees/ -H="2DAlphabet API" -f -e --module-first -d 3 --templatedir source/_templates/ ../TwoDAlphabet
make html