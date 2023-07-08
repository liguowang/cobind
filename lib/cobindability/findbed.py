#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 11:34:34 2021

@author: m102324
"""
import os


def findBedFiles(root_path='.', base_name=False, skip_empty=True):
    """
    Search bed or bed-like files. Supported file formats:
    'bed', 'bed3','bed4','bed6','bed12','bedgraph','bgr','narrowpeak',
    'broadpeak','gappedpeak', 'bb',' bigbed'

    Parameters
    ----------
    root_path : str, optional
        Name of directory to seaerch. The default is the current directory.
    base_name : bool, optional
        Return only the base names of BED files. The default is False.
    skip_empty : bool, optional
        Skip empty BED files. The default is True.

    Returns
    -------
    bed_files : list
        A list of bed or bed-like files.

    """
    suffix_lib = ['bed', 'bed3', 'bed4', 'bed6', 'bed12', 'bedgraph',
                  'bgr', 'narrowpeak', 'broadpeak', 'gappedpeak',
                  'bb', 'bigbed']
    bed_files = []
    for root, dirs, files in os.walk(root_path, followlinks=True):
        for f in files:
            suf = f.split('.')[-1]
            if suf.lower() in suffix_lib:
                long_name = os.path.join(root, f)
                if skip_empty and (os.stat(long_name).st_size == 0):
                    continue
                if base_name:
                    bed_files.append(f)
                else:
                    bed_files.append(long_name)
    return bed_files
