import numpy as np
import matplotlib.pylab as plt
import lmfit
import sys
import matplotlib
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from astropy.stats import sigma_clip
from scipy.interpolate import interp1d
from matplotlib.colors import Normalize
#from numpy.linalg import LinAlgError

from poly import legendre
import fit_params
from fit_params import fit
from phi_bar_sky import pa_bar_sky


from M_tabulated import M_tab


def bisym_mod(vel, evel, guess0, vary, n_it, rstart, rfinal, ring_space, frac_pixel, delta, pixel_scale,bar_min_max,errors, config, e_ISM):
#def bisym_mod(vel, evel, guess0, vary, n_it=5,rstart = 2, rfinal = 55, ring_space = 2, frac_pixel = 0.7, r_back = 20, delta = 1, pixel_scale = 0.2, r_bar_max = 20):

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "bisymmetric"
		[ny,nx] = vel.shape
		shape = [ny,nx]
		r_bar_min, r_bar_max = bar_min_max



		"""

		 					BYSIMETRIC MODEL


		"""

		chisq_global = 1e10
		PA, INC, XC,YC,VSYS = 0,0,0,0,0
		Vrot, Vrad, Vsys,Vtan = [],[],[],[]
		R = 0
		best_xy_pix = []


		rings = np.arange(rstart, rfinal, ring_space)

		for it in np.arange(n_it):
			guess = [vrot0,vr20,pa0,inc0,x0,y0,vsys0,0,theta_b]

			vrot_tab,vrad_tab,vtan_tab = np.asarray([]),np.asarray([]),np.asarray([])
			index = 0

			nrings = len(rings)
			n_annulus = nrings - 1

			
			R_pos = np.asarray([])

			for ring in rings:

				from pixel_params import pixels
				fpix = pixels(shape,vel,pa0,inc0,x0,y0,ring, delta=delta,pixel_scale = pixel_scale)

				if fpix > frac_pixel:


					if ring >=  r_bar_min and ring <= r_bar_max:

							# Create bisymetric model
						try:

							v_rot_k,v_rad_k,v_tan_k = M_tab(pa0,inc0,x0,y0,theta_b,ring, delta,index, shape, vel-vsys0, pixel_scale=pixel_scale,vmode = vmode)

							vrot_tab = np.append(vrot_tab,v_rot_k)
							vrad_tab = np.append(vrad_tab,v_rad_k)
							vtan_tab = np.append(vtan_tab,v_tan_k)
							R_pos = np.append(R_pos,ring)


						except(np.linalg.LinAlgError):
								vrot_tab,vrad_tab,vtan_tab = np.append(vrot_tab,100),np.append(vrad_tab,10),np.append(vtan_tab,10)							
								R_pos = np.append(R_pos,ring)
					else:

							# Create ciruclar model

							v_rot_k,v_rad_k,v_tan_k = M_tab(pa0,inc0,x0,y0,theta_b,ring, delta,index, shape, vel-vsys0, pixel_scale=pixel_scale,vmode = "circular")
							vrot_tab = np.append(vrot_tab,v_rot_k)
							vrad_tab = np.append(vrad_tab,1)
							vtan_tab = np.append(vtan_tab,1)
							R_pos = np.append(R_pos,ring)
	
			guess = [vrot_tab,vrad_tab,pa0,inc0,x0,y0,vsys0, vtan_tab,theta_b]
			v_2D_mdl, vrot , vrad, vsys0,  pa0, inc0, x0, y0, vtan,theta_b, xi_sq, n_data, Errors = fit(shape, vel, evel, guess, vary, vmode, config, R_pos, fit_method = "Powell", e_ISM = e_ISM, pixel_scale = pixel_scale, ring_space = ring_space )


			if xi_sq < chisq_global:

				PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
				Vrot = vrot
				Vrad = vrad
				Vtan = vtan
				chisq_global = xi_sq
				V2D = v_2D_mdl
				Rings = R_pos
				std_errors = Errors
				GUESS = [vrot_tab, vrad_tab, PA, INC, XC, YC, VSYS, vtan_tab, THETA]
				


		# Maybe the code has found the bar minor axis, in such case Vrad <0 and Vtan <0


		if np.nanmean(Vrad) <0 or np.nanmean(Vtan) < 0:
			
			GUESS[-1] = GUESS[-1] - 90
			NEW = GUESS[-1]
			if NEW < 0: GUESS[-1] = 180 + NEW

			VARY = [True, True, False,False,False,False,False,True,True]
			v_2D_mdl_, vrot_ , vrad_,_vsys0_,  pa0_, inc0_, x0_, y0_, vtan_, theta_b_, xi_sq_, n_data_, Errors_ = fit(shape, vel, evel, GUESS, VARY, vmode, config, Rings, fit_method = "Powell", e_ISM = e_ISM, pixel_scale = pixel_scale, ring_space = ring_space )
			if xi_sq_ < chisq_global:
				Vrot = vrot_
				Vrad = vrad_
				Vtan = vtan_
				THETA = theta_b_
				V2D = v_2D_mdl_
				std_errors = Errors_





		Vtan = np.asarray(Vtan)
		Vrad = np.asarray(Vrad)
		Vrot = np.asarray(Vrot)
		PA_bar_major = pa_bar_sky(PA,INC,THETA)
		PA_bar_minor = pa_bar_sky(PA,INC,THETA-90)

		if errors == 1:
			from mcmc import fit_mcmc
			res_mcmc =  fit_mcmc(shape, vel, evel, GUESS, vary, vmode, "", Rings, fit_method = "emcee", e_ISM = e_ISM, pixel_scale = pixel_scale, ring_space = ring_space  )
		else:
			std_Vrot,std_Vrad,std_pa, std_inc, std_x0, std_y0, std_Vsys, std_theta, std_Vtan = std_errors

			res_mcmc = Vrot*0,[std_Vrot,std_Vrot],Vrot*0,[std_Vrad,std_Vrad],0,[std_pa,std_pa],0,[std_inc,std_inc],0,[std_x0,std_x0],0,[std_y0,std_y0],0,[std_Vsys,std_Vsys],0,[std_theta,std_theta],Vrot*0,[std_Vtan,std_Vtan]


		return PA,INC,XC,YC,VSYS,THETA,Rings,Vrot,Vrad,Vtan,V2D,PA_bar_major,PA_bar_minor,chisq_global,res_mcmc
