sudo add-apt-repository -y ppa:cran/poppler
sudo sudo apt-get install -y libpoppler-cpp-dev
mkdir -p ~/.R
echo "CXX14FLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function  -Wno-macro-redefined" >> ~/.R/Makevars
echo "CXX14=g++ -std=c++1y -fext-numeric-literals -fPIC" >> ~/.R/Makevars
Rscript -e 'install.packages(c("devtools", "ndjson"))'
