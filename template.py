import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp1d

M_I = np.genfromtxt("/home/carlos/XS-code/template_RC_optical_radius.csv", usecols = 0)
V0 = np.genfromtxt("/home/carlos/XS-code/template_RC_optical_radius.csv", usecols = 1)
r_RE_r_opt = np.genfromtxt("/home/carlos/XS-code/template_RC_optical_radius.csv", usecols = 2)
alpha = np.genfromtxt("/home/carlos/XS-code/template_RC_optical_radius.csv", usecols = 3)

"""
M_I = np.genfromtxt("template_RC_disk_scale.csv", usecols = 0)
V0 = np.genfromtxt("template_RC_disk_scale.csv", usecols = 1)
r_RE_r_opt = np.genfromtxt("template_RC_disk_scale.csv", usecols = 2)
alpha = np.genfromtxt("template_RC_disk_scale.csv", usecols = 3)
"""



def params(M_I_new):

	f0 = interp1d(M_I, V0)
	f1 = interp1d(M_I, r_RE_r_opt)
	f2 = interp1d(M_I, alpha)
	
	V0_new = f0(M_I_new)
	r_new = f1(M_I_new)
	alpha_new = f2(M_I_new)

	#print(V0_new,r_new,alpha_new)
	return V0_new,r_new,alpha_new


"""
plt.plot(M_I, alpha)
plt.plot(-22.7,params(-22.7)[-1],"ko")
plt.ylim(-0.05,0.15)
plt.gca().invert_xaxis()
plt.show()
"""
