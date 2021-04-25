import numpy as np
import numpy as np
import matplotlib.pylab as plt
import scipy
import sys
import lmfit
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
from matplotlib.gridspec import GridSpec



from weights_interp import weigths_w


 
def Rings(xy_mesh,pa,inc,x0,y0):
	(x,y) = xy_mesh

	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))



	R= np.sqrt(X**2+(Y/np.cos(inc))**2)

	return R


def RADIAL_MODEL2(w1,Vsys,w2,Vrot,w3,Vrad):
	vlos = w1*Vsys + w2*Vrot + w3*Vrad
	return np.ravel(vlos)



def polynomial(x,a0,a1,a2,a3,a4,a5):
	x = np.asarray(x)
	y = a0 + a1*x +a2*x**2 +a3*x**3 + a4*x**4 + a5*x**5
	return y

def linear(x,a0,a1):
	x = np.asarray(x)
	y = a0 + a1*x 
	return y

def fit_polynomial(x,dato):
	x = np.asarray(x)
	dato = np.asarray(dato)


	def residual_line(pars,x,data=None):
		parvals = pars.valuesdict()
		a0 = parvals['a0']
		a1 = parvals['a1']
		a2 = parvals['a2']
		a3 = parvals['a3']
		a4 = parvals['a4']
		a5 = parvals['a5']
		model = polynomial(x,a0,a1,a2,a3,a4,a5)
		objective = model - data
		return objective**2


	fit_param = Parameters()
	fit_param.add('a0', value=0)
	fit_param.add('a1', value=0)
	fit_param.add('a2', value=0)
	fit_param.add('a3', value=0)
	fit_param.add('a4', value=0)
	fit_param.add('a5', value=0)

	out = minimize(residual_line, fit_param, args=(x,), kws={'data': dato},method='Powell', nan_policy = "omit")#, method = "emcee")
	best = out.params
	a0, a1, a2, a3, a4, a5 = best["a0"].value,best["a1"].value, best["a2"].value, best["a3"].value, best["a4"].value, best["a5"].value
	return a0, a1, a2, a3, a4, a5


def fit_linear(x,dato):
	x = np.asarray(x)
	dato = np.asarray(dato)


	def residual_line(pars,x,data=None):
		parvals = pars.valuesdict()
		a0 = parvals['a0']
		a1 = parvals['a1']
		model = linear(x,a0,a1)
		objective = model - data
		return objective**2


	fit_param = Parameters()
	fit_param.add('a0', value=1e3)
	fit_param.add('a1', value=0)

	out = minimize(residual_line, fit_param, args=(x,), kws={'data': dato},method='Powell')#, method = "emcee")
	best = out.params
	a0, a1 = best["a0"].value,best["a1"].value
	return a0, a1


def Errors(res,p,vmode):
	Vrot,Vrad,Vtan = [],[],[]
	eVrot_l,eVrot_u,eVrad_l,eVrad_u,eVtan_l,eVtan_u = [],[],[],[],[],[]


	for name in p.keys():
			
				if "Vsys" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					Vsys = median
					one_sigma_l_Vsys, one_sigma_u_Vsys = abs(err_m1), abs(err_p1)
				if "pa" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					pa = median
					one_sigma_l_pa, one_sigma_u_pa = abs(err_m1), abs(err_p1)

				if "inc" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					inc = median
					one_sigma_l_inc, one_sigma_u_inc = abs(err_m1), abs(err_p1)


				if "x0" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					x0 = median
					one_sigma_l_x0, one_sigma_u_x0 = abs(err_m1), abs(err_p1)


				if "y0" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					y0 = median
					one_sigma_l_y0, one_sigma_u_y0 = abs(err_m1), abs(err_p1)

				if "phi_b" in name:
					quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
					median = quantiles[2]
					err_m2 = quantiles[0] - median # two_sigma_l
					err_m1 = quantiles[1] - median # one_sigma_l
					err_p1 = quantiles[3] - median # one_sigma_u
					err_p2 = quantiles[4] - median # two_sigma_u

					theta = median
					etheta_l, etheta_u = abs(err_m1), abs(err_p1)
		


				if "Vrot" in name:

						quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
						median = quantiles[2]
						err_m2 = quantiles[0] - median # two_sigma_l
						err_m1 = quantiles[1] - median # one_sigma_l
						err_p1 = quantiles[3] - median # one_sigma_u
						err_p2 = quantiles[4] - median # two_sigma_u

						vrot = median
						one_sigma_l_vrot, one_sigma_u_vrot = abs(err_m1), abs(err_p1)
						Vrot.append(vrot)
						eVrot_l.append(one_sigma_l_vrot)
						eVrot_u.append(one_sigma_u_vrot)

				if "Vrad" in name:

						quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
						median = quantiles[2]
						err_m2 = quantiles[0] - median # two_sigma_l
						err_m1 = quantiles[1] - median # one_sigma_l
						err_p1 = quantiles[3] - median # one_sigma_u
						err_p2 = quantiles[4] - median # two_sigma_u

						vrad = median
						one_sigma_l_vrad, one_sigma_u_vrad = abs(err_m1), abs(err_p1)
						Vrad.append(vrad)
						eVrad_l.append(one_sigma_l_vrad)
						eVrad_u.append(one_sigma_u_vrad)

				if "Vtan" in name:

						quantiles = np.percentile(res.flatchain[name], [2.275, 15.865, 50, 84.135, 97.275])
						median = quantiles[2]
						err_m2 = quantiles[0] - median # two_sigma_l
						err_m1 = quantiles[1] - median # one_sigma_l
						err_p1 = quantiles[3] - median # one_sigma_u
						err_p2 = quantiles[4] - median # two_sigma_u

						vtan = median
						one_sigma_l_vtan, one_sigma_u_vtan = abs(err_m1), abs(err_p1)
						Vtan.append(vtan)
						eVtan_l.append(one_sigma_l_vtan)
						eVtan_u.append(one_sigma_u_vtan)

		
	if vmode == "circular": 
		theta, etheta_l,etheta_u = 0,0,0
		eVtan_l,eVtan_u = [],[]
		eVrad_l,eVrad_u = [],[]
	if vmode == "radial": 
		theta, etheta_l,etheta_u = 0,0,0
		eVtan_l,eVtan_u = [],[]



	return Vrot,[eVrot_l,eVrot_u], Vrad,[eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0],Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys],theta,[etheta_l,etheta_u], Vtan,[eVtan_l,eVtan_u]




from kin_components import CIRC_MODEL
from kin_components import RADIAL_MODEL
from kin_components import BISYM_MODEL




def fit_mcmc(shape, vel_map, e_vel_map, guess,vary,vmode,config, rings_pos, ring_space, fit_method = "leastsq", e_ISM = 5, pixel_scale = 1):
	vary = [True,True,True,True,True,True,True,True,True]

	vrot0,vr20,pa0,inc0,X0,Y0,vsys0,vt20,theta0 = guess
	constant_params = [pa0,inc0,X0,Y0,vsys0,theta0]

	[ny,nx] = shape

	nrings = len(rings_pos)
	n_annulus = nrings - 1  


	interp_model = np.zeros((ny,nx))
	mask_ring = np.ones((ny,nx), dtype = bool)



	X = np.arange(0, nx, 1)
	Y = np.arange(0, ny, 1)
	XY_mesh = np.meshgrid(X,Y)

	r_n = Rings(XY_mesh,pa0*np.pi/180,inc0*np.pi/180,X0,Y0)




	if vmode == "circular":


		def vmodel_dataset(pars, i,  xy_mesh,  r_0 = None, r_space = ring_space):


			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% i]
			x0,y0 = parvals['x0'],parvals['y0']

			modl = (CIRC_MODEL(xy_mesh,Vrot,pa,inc,x0,y0))*weigths_w(xy_mesh,shape,pa,inc,x0,y0,r_0,r_space,pixel_scale=pixel_scale)


			return modl,Vsys

		def residual(pars, data= None):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			resid = np.array([])


			# make residual per data set

			for N in range(n_annulus):

				mdl_ev = 0
				r_space_k = rings_pos[N+1] - rings_pos[N] 
				mask = np.where( (r_n >= rings_pos[N] ) & (r_n < rings_pos[N+1]) )
				x,y = XY_mesh[0][mask], XY_mesh[1][mask] 
				XY = (x,y)




				for kk in range(2):
					Vxy,Vsys = vmodel_dataset(pars, N+kk,  XY, r_0 = rings_pos[N], r_space = r_space_k )

					mdl_ev = mdl_ev + Vxy[kk]


					if N == 0 and kk == 0:
						
						mask1 = np.where( (r_n < rings_pos[0] ) )
						x1,y1 = XY_mesh[0][mask1], XY_mesh[1][mask1] 
						XY1 = (x1,y1)


						#
						#
						# inner interpolation
						#
						#
					
						#(a) velocity rise linearly from zero

						r_space_0 =  rings_pos[0]
						Vxy,Vsys = vmodel_dataset(pars, 0,  XY1, r_0 = 0, r_space = r_space_0)
					
						interp_model[mask1] = Vxy[1] + Vsys

				interp_model[mask] = mdl_ev + Vsys




			
			sigma = np.sqrt(e_vel_map**2 + e_ISM**2)
			interp_model[interp_model == 0]  = np.nan


			convolved = 0

			if convolved == True:

				interp_model_conv = convolve(interp_model, gauss_kern, nan_treatment = "interpolate", boundary = "extend", preserve_nan = 1)
			else:
				interp_model_conv = interp_model


			Resid = (vel_map - interp_model_conv)/sigma


			return Resid.flatten()

			



		fit_params = Parameters()


		for iy  in range(nrings):
				fit_params.add('Vrot_%i' % (iy),value=vrot0[iy], vary = vary[0], min = -200, max = 400)

		if config == "":

			fit_params.add('Vsys', value=vsys0, vary = vary[6])
			fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
			fit_params.add('inc', value=inc0, vary = vary[3], min = -90, max = 90)
			fit_params.add('x0', value=X0, vary = vary[4],  min = 0, max = nx)
			fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = ny)


		else:
			k = 0
			for res,const in zip(config_file(config),constant_params):
				#param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), float(res["val"]), float(res["min"]), float(res["max"]) 
				param, fit, val, vmin, vmax = str(res["param"]), bool(float(res["fit"])), eval(res["val"]), eval(res["min"]), eval(res["max"]) 
				fit_params.add(param, value = constant_params[k], vary = fit, min = vmin, max = vmax)
				k = k+1



		nvarys = 2*len(vrot0) + 5
		walkers =5*nvarys
		res = minimize(residual, fit_params,method = fit_method, burn = 300, steps = 1000, thin = 100, nwalkers = walkers,progress = False,nan_policy = "omit")#, nan_policy = 

		Vrot,[eVrot_l,eVrot_u], Vrad,[eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l, etheta_u], Vtan,[eVtan_l,eVtan_u]  = Errors(res,fit_params,vmode)



		return Vrot,[eVrot_l,eVrot_u], Vrad, [eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l,etheta_l], Vtan,[eVtan_l,eVtan_u] 





	if vmode == "radial":
		def vmodel_dataset(pars, i,  xy_mesh,  r_0 = None, r_space = ring_space):


			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% i]
			Vrad = parvals['Vrad_%i'% i]
			x0,y0 = parvals['x0'],parvals['y0']

			modl = (RADIAL_MODEL(xy_mesh,Vrot,Vrad,pa,inc,x0,y0))*weigths_w(xy_mesh,shape,pa,inc,x0,y0,r_0,r_space,pixel_scale=pixel_scale)


			return modl,Vsys

		def residual(pars, data= None):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			resid = np.array([])


			# make residual per data set

			for N in range(n_annulus):

				mdl_ev = 0
				r_space_k = rings_pos[N+1] - rings_pos[N] 
				mask = np.where( (r_n >= rings_pos[N] ) & (r_n < rings_pos[N+1]) )
				x,y = XY_mesh[0][mask], XY_mesh[1][mask] 
				XY = (x,y)




				for kk in range(2):
					Vxy,Vsys = vmodel_dataset(pars, N+kk,  XY, r_0 = rings_pos[N], r_space = r_space_k)

					mdl_ev = mdl_ev + Vxy[kk]


					if N == 0 and kk == 0:
						
						mask1 = np.where( (r_n < rings_pos[0] ) )
						x1,y1 = XY_mesh[0][mask1], XY_mesh[1][mask1] 
						XY1 = (x1,y1)


						#
						#
						# inner interpolation
						#
						#
					
						#(a) velocity rise linearly from zero

						r_space_0 =  rings_pos[0]
						Vxy,Vsys = vmodel_dataset(pars, 0,  XY1, r_0 = 0, r_space = r_space_0)
					
						interp_model[mask1] = Vxy[1] + Vsys

				interp_model[mask] = mdl_ev + Vsys




			
			sigma = np.sqrt(e_vel_map**2 + e_ISM**2)
			interp_model[interp_model == 0]  = np.nan


			convolved = 0

			if convolved == True:

				interp_model_conv = convolve(interp_model, gauss_kern, nan_treatment = "interpolate", boundary = "extend", preserve_nan = 1)
			else:
				interp_model_conv = interp_model


			Resid = (vel_map - interp_model_conv)/sigma


			return Resid.flatten()


    
#####




		fit_params = Parameters()

		for iy  in range(nrings):

				fit_params.add('Vrot_%i' % (iy),value=vrot0[iy], vary = vary[0], min = 0, max = 400)
				fit_params.add('Vrad_%i' % (iy), value=vr20[iy], vary = vary[1], min = -200, max = 200)



		fit_params.add('Vsys', value=vsys0, vary = vary[6])
		fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
		fit_params.add('inc', value=inc0, vary = vary[3], min = 0, max = 90)
		fit_params.add('x0', value=X0, vary = vary[4], min = 0, max = nx)
		fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = ny)



		nvarys = 2*len(vrot0) + 5
		walkers =5*nvarys

		res = minimize(residual, fit_params,method = fit_method, burn = 300, steps = 2000, thin = 100, nwalkers = walkers, progress = False, nan_policy = "omit") 

		Vrot,[eVrot_l,eVrot_u], Vrad,[eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l, etheta_u], Vtan,[eVtan_l,eVtan_u]  = Errors(res,fit_params,vmode)


		return Vrot,[eVrot_l,eVrot_u], Vrad, [eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l,etheta_l], Vtan,[eVtan_l,eVtan_u] 









	if vmode == "bisymmetric":

		def vmodel_dataset(pars, i,  xy_mesh,  r_0 = None, r_space = ring_space):


			parvals = pars.valuesdict()
			pa = parvals['pa']
			inc = parvals['inc']
			Vsys = parvals['Vsys']
			Vrot = parvals['Vrot_%i'% i]
			Vrad = parvals['Vrad_%i'% i]
			Vtan = parvals['Vtan_%i'% i]
			x0,y0 = parvals['x0'],parvals['y0']
			phi_b = parvals['phi_b']

			modl = (BISYM_MODEL(xy_mesh,Vrot,Vrad,pa,inc,x0,y0,Vtan,phi_b))*weigths_w(xy_mesh,shape,pa,inc,x0,y0,r_0,r_space,pixel_scale=pixel_scale)


			return modl,Vsys

		def residual(pars, data= None):
			"""Calculate total residual for fits of VMODELS to several data sets."""

			resid = np.array([])


			# make residual per data set

			for N in range(n_annulus):

				mdl_ev = 0
				r_space_k = rings_pos[N+1] - rings_pos[N] 
				mask = np.where( (r_n >= rings_pos[N] ) & (r_n < rings_pos[N+1]) )
				x,y = XY_mesh[0][mask], XY_mesh[1][mask] 
				XY = (x,y)




				for kk in range(2):
					Vxy,Vsys = vmodel_dataset(pars, N+kk,  XY, r_0 = rings_pos[N], r_space = r_space_k)

					mdl_ev = mdl_ev + Vxy[kk]


					if N == 0 and kk == 0:
						
						mask1 = np.where( (r_n < rings_pos[0] ) )
						x1,y1 = XY_mesh[0][mask1], XY_mesh[1][mask1] 
						XY1 = (x1,y1)


						#
						#
						# inner interpolation
						#
						#
					
						#(a) velocity rise linearly from zero

						r_space_0 =  rings_pos[0]
						Vxy,Vsys = vmodel_dataset(pars, 0,  XY1, r_0 = 0, r_space = r_space_0)
					
						interp_model[mask1] = Vxy[1] + Vsys

				interp_model[mask] = mdl_ev + Vsys




			
			sigma = np.sqrt(e_vel_map**2 + e_ISM**2)
			interp_model[interp_model == 0]  = np.nan


			convolved = 0

			if convolved == True:

				interp_model_conv = convolve(interp_model, gauss_kern, nan_treatment = "interpolate", boundary = "extend", preserve_nan = 1)
			else:
				interp_model_conv = interp_model


			Resid = (vel_map - interp_model_conv)/sigma


			return Resid.flatten()




		fit_params = Parameters()


		for iy in range(nrings):
				if vr20[iy]-1 == 0 and vt20[iy]-1 ==0:
					vary_vrad = False
					vary_vtan = False


				else:

					vary_vrad = vary[1]
					vary_vtan = vary[7]
					#vary_theta = vary[8]


				fit_params.add('Vrot_%i' % (iy),value=vrot0[iy], vary = vary[0], min = -200, max = 400)
				fit_params.add('Vrad_%i' % (iy), value=vr20[iy], vary = vary_vrad,  min = -200, max = 200)
				fit_params.add('Vtan_%i' % (iy), value=vt20[iy], vary = vary_vtan, min = -200, max = 200)



		fit_params.add('Vsys', value=vsys0, vary = vary[6])
		fit_params.add('pa', value=pa0, vary = vary[2], min = 0, max = 360)
		fit_params.add('inc', value=inc0, vary = vary[3], min = 0, max = 90)
		fit_params.add('x0', value=X0, vary = vary[4], min = 0, max = nx)
		fit_params.add('y0', value=Y0, vary = vary[5], min = 0, max = nx)
		fit_params.add('phi_b', value=theta0, vary = vary[7], min = 0 , max = 360)



		res = minimize(residual_bisym, fit_params,method = fit_method, burn = 300, steps = 2000, nwalkers = 160,thin = 100, progress = False, nan_policy = "omit")

		Vrot,[eVrot_l,eVrot_u], Vrad,[eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l, etheta_u], Vtan,[eVtan_l,eVtan_u]  = Errors(res,fit_params,vmode)


		return Vrot,[eVrot_l,eVrot_u], Vrad, [eVrad_l,eVrad_u], pa,[one_sigma_l_pa, one_sigma_u_pa], inc,[one_sigma_l_inc, one_sigma_u_inc], x0,[one_sigma_l_x0, one_sigma_u_x0], y0,[one_sigma_l_y0, one_sigma_u_y0], Vsys,[one_sigma_l_Vsys, one_sigma_u_Vsys], theta,[etheta_l,etheta_l], Vtan,[eVtan_l,eVtan_u] 



