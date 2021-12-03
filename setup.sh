#install python libraries
pip3 install -r import_doc/requirements.txt
#install BiocManager and topGO
sudo apt-get install libcurl4-openssl-dev libssl-dev
sudo R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); library('BiocManager'); BiocManager::install('topGO'); BiocManager::install('Rgraphviz')"
