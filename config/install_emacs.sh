#!/bin/bash

mkdir comp
cd comp
wget http://mirror.switch.ch/ftp/mirror/gnu/emacs/emacs-24.5.tar.gz -O -|tar xzf -
cd emacs*
./configure --prefix=$HOME
make -j 2
make install
cd ..
rm -fr emacs*
cd ..
rm -fr comp

echo "Add:
(savehist-mode 1)
(add-to-list 'auto-mode-alist '(\"\\.ypp\\'\" . c++-mode))
(add-to-list 'auto-mode-alist '(\"\\.cu\\'\" . c++-mode))
(add-to-list 'auto-mode-alist '(\"\\.lpp\\'\" . c++-mode))
 to ~/.emacs"
