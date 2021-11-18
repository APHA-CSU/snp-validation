import matplotlib.pyplot as plt
import numpy as np

TP = 1
E = np.linspace(0, 10*TP, 100)

F = TP / (TP + 0.5*E)
A = TP / (TP + E)

plt.plot(E, F, label="F-score")
plt.plot(E, A, label="Accuracy")

plt.xlabel('(FP+FN)')
plt.ylabel('Metric')
plt.ylim([0, 1])
plt.legend()
plt.grid()

plt.show()
