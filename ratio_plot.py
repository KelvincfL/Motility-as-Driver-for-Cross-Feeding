import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

ldhcd = [1,0.85,0.3,0.2,0]
ldhcc = [0,0.1,0.6,0.8,1]
ldmcd = [1,1,1,0.55,0.55]
ldmcc = [0,0,0,0.45,0.45]
ldlcd = [1,1,0.95,0.7,0.55]
ldlcc = [0,0,0,0.3,0.45]

hdhcd = [0,0,0,0,0]
hdhcc = [1,1,1,1,1]
hdmcd = [0.15,0.05,0,0,0]
hdmcc = [0.35,0.95,1,1,1]
hdlcd = [0.4,0.1,0,0,0]
hdlcc = [0.35,0.9,1,1,1]

biased_ldhcd = [0.625,0.625,0.625,0.375]
biased_ldhcc = [0.25,0.25,0.25,0.375]
biased_ldmcd = [1,1,1,1]
biased_ldmcc = [0,0,0,0]
biased_ldlcd = [1,1,1,1]
biased_ldlcc = [0,0,0,0]

antibiased_ldhcd = [0.375,0.5,0.75,1]
antibiased_ldhcc = [0.25,0.25,0.125,0]
antibiased_ldmcd = [1,1,1,1]
antibiased_ldmcc = [0,0,0,0]
antibiased_ldlcd = [1,1,1,1]
antibiased_ldlcc = [0,0,0,0]

biased_range = np.logspace(-5, -2, base=10, num=4, endpoint=True)
ratio_range = [0.1,0.5,1,2,10]
sample_size = 8
D5error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldlcd])
C5error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldlcc])

D10error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldmcd])
C10error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldmcc])

D50error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldhcd])
C50error = np.array([stats.norm.ppf(0.95)*np.sqrt(i*(1-i)/sample_size) for i in antibiased_ldhcc])

plt.figure(figsize=(30, 40))
plt.errorbar(biased_range,antibiased_ldlcd, D5error, linestyle = '-.', label='Defector: Cost 5%',color= 'royalblue',marker='o')
plt.errorbar(biased_range,antibiased_ldlcc, C5error, linestyle = '-.', label='Cooperator: Cost 5%',color= 'gold',marker='o')
plt.errorbar(biased_range,antibiased_ldmcd, D10error, linestyle = '--', label='Defector: Cost 10%',color= 'blue',marker='o')
plt.errorbar(biased_range,antibiased_ldmcc, C10error, linestyle = '--', label='Cooperator: Cost 10%',color= 'orange',marker='o')
plt.errorbar(biased_range,antibiased_ldhcd, D50error,linestyle = '-', label='Defector: Cost 50%',color= 'darkblue',marker='o')
plt.errorbar(biased_range,antibiased_ldhcc, C50error,linestyle = '-',label='Cooperator: Cost 50%',color= 'darkorange',marker='o')
plt.grid(True)
plt.xscale('log')
plt.xlabel('Exchange Rate')
plt.ylabel('Extinction Probability')
plt.title('Biased Motility, Low Density, Chemoaversion')
plt.legend()
plt.show()