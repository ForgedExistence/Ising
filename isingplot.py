import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

particle = 100*100
data_meas = np.loadtxt('equi_test.dat')
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

def exp_decay(x, a, b):
    return a * np.exp(-b * x)

plt.figure()
plt.title("auto-cor")
time = np.arange(0, tmax, step)
plt.grid()
plt.scatter(time, np.log(auto_cor_normalized), s = 1)
plt.show()



y = np.log(auto_cor_normalized)

y[np.isnan(y)] = np.nanmean(y)
y[np.isinf(y)] = np.nanmedian(y)

popt, pcov = curve_fit(exp_decay, time, y)

a_opt = popt[0]
b_opt = popt[1]

x = (-np.log(0.5) / b_opt)
in_x = (1/x)*tmax

print(in_x)

