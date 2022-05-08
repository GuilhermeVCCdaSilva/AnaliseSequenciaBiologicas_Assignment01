sudo apt-get update

#PAUP
wget http://phylosolutions.com/paup-test/paup4a168_ubuntu64.gz
gunzip paup4a168_ubuntu64.gz
chmod 755 paup4a168_ubuntu64
sudo mv paup4a168_ubuntu64 /usr/local/bin/paup
sudo apt -y install libpython2.7

#FigTree
sudo apt -y install figtree

#Mafft
sudo apt -y install mafft

#Seqmagik
sudo apt-get -y install seqmagick

#Possível erro biopython com seqmagik
pip install biopython==1.77

#toytree
pip3 install toytree

#python--dependencies
pip install biopython
pip install biopython --upgrade
sudo apt install python3-pip
pip3 --version
pip3 install biopython