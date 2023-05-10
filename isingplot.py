import numpy as np
import matplotlib.pyplot as plt

particle = 100*100
file = open("equi_test.dat",  "r")
data_meas = np.genfromtxt(file)
file.close()
average_mag = data_meas[:, 1]
print(average_mag[1])
tmax = 59900000
step = 100000


def auto_correlation_function(t):
    factor = 1/(tmax-t*step)
    sum_1 = 0
    sum_2 = 0
    sum_3 = 0
    t_prime = np.arange(0, int((tmax/step-t)))
    for i in t_prime:
        sum_1 += average_mag[i]*average_mag[t+i]
        sum_2 += average_mag[i]
        sum_3 += average_mag[t+i]

    return factor*sum_1-factor*sum_2*factor*sum_3


auto_cor = []
for i in range(0, int(tmax/step), 1):
    auto_cor.append(auto_correlation_function(i))
auto_cor_normalized = np.array(auto_cor)/auto_cor[0]
plt.figure()
plt.title("auto-cor")
time = np.arange(0, tmax, step)
plt.grid()
plt.plot(time, np.log(auto_cor_normalized), '.')
plt.show()
