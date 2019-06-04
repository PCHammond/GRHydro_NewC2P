#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:51:16 2019

@author: peter
"""

for i in range(0,48):
    print((i%3)+1, (i//3)%4, (i//12)%4,((i%3)+1) + 4*((i//3)%4) + 16*((i//12)%4) + 1)

for i in range(0,48):
    print((i%4), (i//4)%3 + 1, (i//12)%4,((i%4)) + 4*((i//4)%3 + 1) + 16*((i//12)%4) + 1)
    
for i in range(0,48):
    print((i%4), (i//4)%4, (i//16)%3 + 1,((i%4)) + 4*((i//4)%4) + 16*((i//16)%3 + 1) + 1)