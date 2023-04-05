#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:17:05 2023

@author: broszms
"""
import sys
from topology import topWriter
from topology import itpGenerator
from topology import mapGoModel
#
def main(file):
    itpGenerator(file)
    mapGoModel()
    topWriter()


if __name__ == '__main__':
    print('TODO')