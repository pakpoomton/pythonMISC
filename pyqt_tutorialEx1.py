#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
ZetCode PyQt4 tutorial 

In this example, we create a simple
window in PyQt4.

author: Jan Bodnar
website: zetcode.com 
last edited: October 2011
"""

import sys
from PyQt4 import QtGui


def main():
    
    app = QtGui.QApplication(sys.argv)

    w = QtGui.QWidget()# intiate GUI window object
    w.resize(250, 150) #size of GUI window
    w.move(300, 300)   #location of GUI window
    w.setWindowTitle('Simple')
    w.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
