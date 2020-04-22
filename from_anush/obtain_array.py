import sys
import csv
import os
import os.path
import pandas as pd
import numpy as np


#program that creates array of genes and data associated with each gene (raw count, length, RPKM) per
#patient seen 

patient = (os.environ.get('file'))
counter = (os.environ.get('p'))



#temp_data = np.loadtxt(patient,skiprows = 1, usecols = (1,2,3))
#print(temp_data.shape)
#temp_data = np.full(shape = temp_data.shape, fill_value=-7)

#look into each patient file and print each array line containing gene, raw count, length, and RPKM
with open(patient) as tab:
    for line in csv.reader(tab, dialect="excel-tab"):
        print(line)
        
