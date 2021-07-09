#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 12:59:37 2021

@author: fabian
"""

import numpy as np
import pandas as pd


#data = np.loadtxt("./mtz.txt")
with open("./mcf.txt") as f:
    content = f.readlines()
content = [x.strip() for x in content] 
    
data = np.array(content).reshape(-1,9)
df = pd.DataFrame(data)

df[[1,2,4,5,8]] = df[[1,2,4,5,8]].astype(str).astype(int)
df[[0,3]] = df[[0,3]].astype(str)
df[[6,7]] = df[[6,7]].astype(str).astype(float)

col_names = ['model_type', 'k', 'total_nodes', 'solution_status',
              'b&b_nodes','best_obj','cpu_time [s]','optimality_gap [%]','b&b_nodes1']

df.columns = col_names

cpu_times = np.array(df[col_names[6]])
cpu_times[1:] = cpu_times[1:]-cpu_times[:-1]
df[col_names[6]] = np.round(cpu_times,0).astype("int")

df[col_names[7]] = np.round(np.array(df[col_names[7]])*100,2)
df[col_names[2]] = np.array(df[col_names[2]])-1

df["instance"] = (np.arange(data.shape[0])/2).astype("int")+1

tab_cols = np.array(["instance"]+col_names)[[0,2,3,4,5,6,7,8]]

print(df[tab_cols].to_latex(index=False))


