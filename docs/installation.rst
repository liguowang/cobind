Installation
=============

cobind is written in Python. Python3 (v3.5.x) is required to run all programs in
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

Install to virtual environment
-------------------------------

Please install `pip <https://pypi.org/project/pip/>`_ (Package Installer for Python) first if you do not have it.


Python's **Virtual Environments** allow Python packages to be installed in an isolated location rather than being installed globally. If you would like to install *cobind* into a virtual environment, please follow `these instructions <https://packaging.python.org/en/latest/tutorials/installing-packages/#creating-and-using-virtual-environments>`_. 
Specifically, follow these steps:

.. code-block::

   $python3 -m venv cobind
   $source cobind/bin/activate
   $pip3 install cobind


Install globally
-----------------
Please install `pip <https://pypi.org/project/pip/>`_ (Package Installer for Python) first if you do not have it.

.. code-block::
   
   #check if pip is available. 
   $ pip --version
   pip 23.0.1 from /Users/m102324/miniconda3/lib/python3.10/site-packages/pip (python 3.10)

Then you can install *cobind* using `pip <https://pypi.org/project/pip/>`_ from `PyPI <https://pypi.org/project/cobind/>`_ or `GitHub <https://github.com/liguowang/cobind>`_

.. code-block::
 
   $ pip install cobind
   #or 
   $ pip install git+https://github.com/liguowang/cobind.git

Upgrade
-------
.. code-block::

   $ pip install cobind --upgrade 

Uninstall
---------

.. code-block::

   $ pip uninstall cobind


