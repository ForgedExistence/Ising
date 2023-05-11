import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


particle = 100*100
for m in np.arange(0.1, 10, 0.1):
    locals()['data_{:.1f}'.format(m)] = np.loadtxt("Results/equi_T={:.1f}.dat".format(m))
#data = np.loadtxt('Results/equi_T=.dat')
    locals()['avg_mag_{:.1f}'.format(m)] = locals()['data_{:.1f}'.format(m)][:,1]
    locals()['auto_cor_{:.1f}'.format(m)] = []
                                                                                    
#print(average_mag[1])
tmax = 99900000-50000000
step = 100000

def auto_correlation_function(average_mag):
    sum_1 = 0
    sum_2 = 0
    sum_3 = 0
    for j in range(0, int(tmax/step), 1):
        factor = 1/(tmax-j*step)
        t_prime = np.arange(0, int((tmax/step-j)))
        for i in t_prime:
            sum_1 += average_mag[i]*average_mag[j+i]
            sum_2 += average_mag[i]
            sum_3 += average_mag[j+i]

    return factor*sum_1-factor*sum_2*factor*sum_3

def straight_line(x, a, b):
    return a*x+b

#auto_cor = []

time = np.arange(0, tmax, step)
cor_time = []

for m in np.arange(0.1, 10, 0.1):
    locals()['auto_cor_{:.1f}'.format(m)].append(auto_correlation_function(locals()['avg_mag_{:.1f}'.format(m)]))
    locals()['auto_cor_norm_{:.1f}'.format(m)] = np.log(np.array(locals()['auto_cor_{:.1f}'.format(m)])/locals()['auto_cor_{:.1f}'.format(m)][0])
    locals()['popt_{:.1f}'.format(m)], locals()['pcov_{:.1f}'.format(m)] = curve_fit(straight_line, time, locals()['auto_cor_norm_{:.1f}'.format(m)])
    cor_time.append(locals()['popt_{:.1f}'.format(m)][0])
    
cor_time = np.array(cor_time)

temp = np.arange(0.1, 10, 0.1)

# plt.scatter(temp, cor_time, s = 0.1)
# plt.gca().set_ylim(auto=True)
# plt.show()
print(cor_time)

        
        
    #auto_cor.append(auto_correlation_function(i))
#auto_cor_normalized = np.array(auto_cor)/auto_cor[0]



# plt.figure()
# plt.title("auto-cor")
# plt.grid()
# plt.scatter(time, np.log(auto_cor_normalized), s=1)
# plt.show()
# y = np.log(auto_cor_normalized)
# y[np.isnan(y)] = np.nanmean(y)
# y[np.isinf(y)] = np.nanmedian(y)

# popt, pcov = curve_fit(straight_line, time, y)

# a_opt = popt[0]
# b_opt = popt[1]

#print(a_opt*step)
