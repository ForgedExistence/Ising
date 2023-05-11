import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


particle = 100*100
tmax = 99900000-50000000
step = 100000
cor_time = []
for m in np.arange(0.1, 10, 0.1):
    data = np.loadtxt("Results/equi_T={:.1f}.dat".format(m))
#data = np.loadtxt('Results/equi_T=.dat')
    av_mag = data[:, 1]
    auto_cor = []
    #print(av_mag)
    def auto_correlation_function(t):
        factor = 1/(tmax-t*step)
        sum_1 = 0
        sum_2 = 0
        sum_3 = 0
        t_prime = np.arange(0, int((tmax/step-t)))
        for i in t_prime:
            sum_1 += av_mag[i]*av_mag[t+i]
            sum_2 += av_mag[i]
            sum_3 += av_mag[t+i]

        return factor*sum_1-factor*sum_2*factor*sum_3
    for i in range(0, int(tmax/step), 1):
        auto_cor.append(auto_correlation_function(i))
    #print(auto_cor)

    auto_cor_norm = np.log(np.array(auto_cor)/auto_cor[0])
    auto_cor_norm[np.isnan(auto_cor_norm)] = 0
    print(auto_cor_norm)
    print(m)
    def straight_line(x, a, b):
        return a*x+b

    time = np.arange(0, tmax, step)

    popt, pcov = curve_fit(straight_line, time, auto_cor_norm)
    cor_time.append(popt[0])

temp = np.arange(0.1, 10, 0.1)

print(cor_time)
plt.figure()
plt.title("Time-cor")
plt.grid()
plt.scatter(temp, cor_time)
plt.show()




#auto_cor = []

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
