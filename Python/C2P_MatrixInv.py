#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:49:40 2019

@author: peter
"""

import numpy as np

def C2P_MatrixInverse(Matrix):
    a = Matrix[0,0]
    b = Matrix[0,1]
    c = Matrix[0,2]
    d = Matrix[1,0]
    e = Matrix[1,1]
    f = Matrix[1,2]
    g = Matrix[2,0]
    h = Matrix[2,1]
    i = Matrix[2,2]
    
    A =  (e*i - f*h)
    B = -(d*i - f*g)
    C =  (d*h - e*g)
    D = -(b*i - c*h)
    E =  (a*i - c*g)
    F = -(a*h - b*g)
    G =  (b*f - c*e)
    H = -(a*f - c*d)
    I =  (a*e - b*d)
    
    det = a*A + b*B + c*C
    
    Inverse = np.array([[A,D,G],
                        [B,E,H],
                        [C,F,I]])/det
    
    return Inverse