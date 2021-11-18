import matplotlib
matplotlib.use('TkAgg')


import matplotlib.pyplot as plt
import pandas as pd


csv_path = '/home/aaronfishman/temp-results/fold-coverage-test-5/stats.csv'

df = pd.read_csv(csv_path)
fold = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50]

ax = plt.figure()

plt.plot(fold, df["FN"])
plt.plot(fold, df["FN"], 'rx')
plt.xlabel('Fold Coverage')
plt.ylabel('False Negatives')
plt.grid()

plt.show()


