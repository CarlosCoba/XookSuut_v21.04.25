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
from M_tabulated import M_tab




def circ_mod(vel, evel, guess0, vary, n_it, rstart, rfinal, ring_space, frac_pixel, delta, pixel_scale, bar_min_max, errors, config, e_ISM ):

		vrot0,vr20,pa0,inc0,x0,y0,vsys0,vtan,theta_b = guess0
		vmode = "circular"
		[ny,nx] = vel.shape
		shape = [ny,nx]




		"""

		 					CIRCULAR MODEL


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

					v_rot_k,v_rad_k,v_tan_k = M_tab(pa0,inc0,x0,y0,theta_b,ring, delta,index, shape, vel-vsys0, pixel_scale=pixel_scale,vmode = vmode)

					vrot_tab = np.append(vrot_tab,v_rot_k)
					R_pos = np.append(R_pos,ring)



			guess = [vrot_tab+1,vrad_tab+1,pa0,inc0,x0,y0,vsys0, vtan_tab,theta_b]
			v_2D_mdl, vrot , vsys0,  pa0, inc0, x0, y0, xi_sq, n_data, Errors = fit(shape, vel, evel, guess, vary, vmode, config, R_pos, fit_method = "powell", e_ISM = e_ISM, pixel_scale = pixel_scale, ring_space = ring_space  )



			
			if xi_sq < chisq_global:

							PA, INC, XC,YC,VSYS,THETA = pa0, inc0, x0, y0,vsys0,theta_b
							Vrot = vrot
							#XY_PIX_POS = xy_pos
							chisq_global = xi_sq
							V2D = v_2D_mdl
							#Rings = rings
							Rings = R_pos
							std_errors = Errors
							GUESS = [Vrot, 0, PA, INC, XC, YC, VSYS, 0, 0]

						



		Vrot = np.array(Vrot)


		if errors == 1:

			from mcmc import fit_mcmc
			res_mcmc =  fit_mcmc(shape, vel, evel, GUESS, vary, vmode, "", Rings, fit_method = "emcee", e_ISM = e_ISM, pixel_scale = pixel_scale, ring_space = ring_space  )

		else:

			std_Vrot,std_Vrad,std_pa, std_inc, std_x0, std_y0, std_Vsys, std_theta, std_Vtan = std_errors

			res_mcmc = Vrot*0,[std_Vrot,std_Vrot],Vrot*0,[std_Vrad,std_Vrad],0,[std_pa,std_pa],0,[std_inc,std_inc],0,[std_x0,std_x0],0,[std_y0,std_y0],0,[std_Vsys,std_Vsys],0,[std_theta,std_theta],Vrot*0,[std_Vtan,std_Vtan]

		return PA,INC,XC,YC,VSYS,0,Rings,Vrot,0*Vrot,0*Vrot,V2D,chisq_global,res_mcmc



