# Install python3
sudo apt-get install -y python3;
# Download pip for **python 3** using wget command:
wget https://bootstrap.pypa.io/get-pip.py -O get-pip.py;
# Install pip for python 3:
sudo python3 get-pip.py;
# Install numpy
python3 -m pip install numpy;
# Install scipy
python3 -m pip install scipy;
# Install matplotlib
sudo apt-get install -y python3-matplotlib;
# Install LaTeX
sudo apt-get install -y texlive;
# Install extra packages for TeX output with matplotlib
sudo apt-get install -y texlive-fonts-recommended texlive-fonts-extra;
sudo apt-get install -y dvipng;
# Install okular, an open-source pdf document reader
sudo apt-get install -y okular;