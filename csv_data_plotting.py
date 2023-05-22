import matplotlib.pyplot as plt
import numpy as np
import csv

f1 = open('scan8863_H2_RHF.csv', 'r')
raw_data1 = [row for row in csv.reader(f1)]

axes_labels1 = raw_data1[0]
data1 = [[float(x), float(y.split(';')[0])] for x,y in raw_data1[1:]]

x_axis, y_axis = np.transpose(data1)
#print(x_axis)
plt.plot(x_axis, y_axis)
plt.xlabel(axes_labels1[0])
plt.ylabel(axes_labels1[1])
#plt.xlim(-1, 0)
#plt.ylim(0.5, 4)

plt.show()

f1.close()