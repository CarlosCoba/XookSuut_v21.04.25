import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from template import params
from scipy import special
from scipy.special import gamma

def Polyex(r, M_I):
	#r_PE = r_PE / r_opt
	#r = r/r_opt
	V0,r_PE,alpha  = params(M_I)

	#r = r/40
	#r_PE = r_PE*40
	V = V0*(1-np.exp(-r/r_PE))*(1+alpha*r/r_PE)
	return V


def xi_sq_pdf(x,k=3):
	xi = (1/2.)**(k/2.)*(np.exp(-x/2.))*x**(k/2.-1)/gamma(k/2.)
	xi = xi/np.nanmax(xi)
	return xi






def plot_true(name,vmode,axis):

	galaxy = "test.%s.fits"%name
	hdu = fits.open("/media/carlos/COBA/XS-code/%s"%galaxy)
	hdr = hdu[0].header
	PA, INC,XC,YC,VSYS,PHI_BAR,VMAX,AMP_VT2,AMP_VR2,M_I,R_OPT,NU1,NU2 = hdr["PA"], hdr["INC"],hdr["XC"],hdr["YC"],hdr["VSYS"],hdr["PHI_BAR"],hdr["VMAX"],hdr["AMP_VT2"],hdr["AMP_VR2"],hdr["M_I"],hdr["R_OPT"],hdr["NU1"],hdr["NU2"] 
	

	R_norm = np.linspace(0,1.8,1000)
	Vrot = Polyex(R_norm, M_I)
	Vt2 = xi_sq_pdf(R_norm*R_OPT,NU1)*AMP_VT2
	Vr2 = xi_sq_pdf(R_norm*R_OPT,NU2)*AMP_VR2
	R1 = R_norm*R_OPT
	R1[R1>R_OPT] = np.nan

	#R = np.linspace(0,38,200)
	#Vrot = 80*(R/27.8)**0.35
	#R1 = R[R < 10*np.pi]
	#R1[R1 == 0] = np.nan
	#Vt2 = 25*np.sin(R1/10.)
	#Vr2 = 0.85*Vt2

	if vmode == "bisymmetric": 

		axis.plot(R1,Vt2,color = "dodgerblue",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R1,Vr2,color = "orange",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R1,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

	if vmode == "radial": 

		axis.plot(R1,Vr2,color = "orange",linestyle='-', alpha = 0.5, linewidth=2)
		axis.plot(R1,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

	if vmode == "circular":
		axis.plot(R1,Vrot,color = "k",linestyle='-', alpha = 0.5, linewidth=2)

