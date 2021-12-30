Installation
=============

cobind are written in Python. Python3 (v3.5.x) is required to run all programs in
cobind.

Dependencies
------------
- `pandas <https://pandas.pydata.org/>`_
- `numpy <http://www.numpy.org/>`_
- `scipy <https://www.scipy.org/>`_
- `bx-python <https://github.com/bxlab/bx-python>`_
- `pyBigWig <https://pypi.org/project/pyBigWig/>`_

.. note::
   These packages will be automatically installed when you use `pip3 <https://pip.pypa.io/en/stable/installing/>`_ to install cobind.

Virtual environment
-------------------
If you would like to install *cobind* into a virtual environment, please follow `these instructions <https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-and-using-virtual-environments>`_. 

Install
-------
Please install `pip <https://pypi.org/project/pip/>`_ (Package Installer for Python) first if you do not have it.

::
 
 #check if pip is avaiable. 
 $ pip --version
 pip 20.1.1 from /Users/m102324/opt/anaconda3/lib/python3.8/site-packages/pip (python 3.8)

Then you can install *cobind* using `pip <https://pypi.org/project/pip/>`_ from `PyPI <https://pypi.org/project/cobind/>`_ or `github <https://github.com/liguowang/cobind>`_

::
 
 $ pip install cobind
 or 
 $ pip install git+https://github.com/liguowang/cobind.git

Upgrade
-------
::

 $ pip install cobind --upgrade 

Uninstall
---------
::

 $ pip uninstall cobind


