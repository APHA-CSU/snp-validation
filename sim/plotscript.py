import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import glob
import os
from pathlib import Path
import shutil


dirs = glob.glob("/home/cameronnicholls/boxplots/OFAT/*/stats.csv")

values = pd.DataFrame(index =[1,2,3,4,5,6,7,8,9,10])
branches=['branch:']

for directory in dirs:
    path = Path(directory)
    sample_name = path.parts[-2]
    branches.append(sample_name)
    csv = pd.read_csv(directory)
    values[sample_name] = csv['FN']

print(branches)
plt.boxplot(values, showfliers = True)
# y = (0.8, 0.85, 0.9, 0.95, 1)
#y = (0.9, 0.92, 0.94, 0.96, 0.98, 1)
plt.xticks(range(0, len(branches)),  branches, rotation=-20)
#plt.yticks(y)
plt.show()
