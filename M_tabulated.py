import numpy as np
import matplotlib.pylab as plt
 
def Rings(xy_mesh,pa,inc,x0,y0,pixel_scale=1):
	(x,y) = xy_mesh


	X = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))
	Y = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))

	R= np.sqrt(X**2+(Y/np.cos(inc))**2)
	R = R*pixel_scale


	return R





def Vlos_BISYM(xy_mesh,Vrot,Vr2,pa,inc,x0,y0,Vsys):
	(x,y) = xy_mesh
	R  = Rings(xy_mesh,pa,inc,x0,y0)
	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)
	vlos = Vsys+np.sin(inc)*(Vrot*cos_tetha + Vr2*sin_tetha)
	return np.ravel(vlos)


def Phi_bar_sky(pa,inc,phi_b_gal):

	pa, inc , phi_b_gal = pa*np.pi/180, inc*np.pi/180 , phi_b_gal*np.pi/180
	phi_sky = pa + np.arctan(np.tan(phi_b_gal)*np.cos(inc))
	return phi_sky*180/np.pi




def weigths_w(xy_mesh,shape,pa,inc,x0,y0,ring,delta,k,pixel_scale):

	r_n = Rings(xy_mesh,pa,inc,x0,y0,pixel_scale)
	a_k = ring


	mask = np.where( (r_n >= a_k - delta) & (r_n < a_k + delta) ) 
	r_n = r_n[mask]

	w_k_n = (1 - (r_n -a_k)/delta)
	w_k_plus_1_n = (r_n -a_k)/delta


	(x,y) = xy_mesh


	return np.ravel(w_k_n),np.ravel(w_k_plus_1_n),mask

def cos_sin(xy_mesh,pa,inc,x0,y0,pixel_scale=1):
	(x,y) = xy_mesh
	R  = Rings(xy_mesh,pa,inc,x0,y0,pixel_scale)

	cos_tetha = (- (x-x0)*np.sin(pa) + (y-y0)*np.cos(pa))/R
	sin_tetha = (- (x-x0)*np.cos(pa) - (y-y0)*np.sin(pa))/(np.cos(inc)*R)


	#return np.ravel(cos_tetha),np.ravel(sin_tetha)
	return cos_tetha,sin_tetha


def trigonometric_weights(xy_mesh,pa,inc,x0,y0,phi_b,mask,vmode="radial",pixel_scale=1):

	cos,sin = cos_sin(xy_mesh,pa,inc,x0,y0)
	cos,sin = cos[mask],sin[mask]

	if vmode == "circular":
		w_rot = np.sin(inc)*cos
		return np.ravel(w_rot)




	if vmode == "radial":
		w_rot = np.sin(inc)*cos
		w_rad = np.sin(inc)*sin
		return np.ravel(w_rot), np.ravel(w_rad)

	if vmode == "bisymmetric":



		theta = np.arctan(sin/cos)
		phi_b = phi_b
		theta_b = theta - phi_b


		w_rot = np.sin(inc)*cos
		w_rad = np.sin(inc)*sin*np.sin(2*theta_b)
		w_tan = np.sin(inc)*cos*np.cos(2*theta_b)
		return np.ravel(w_rot), np.ravel(w_rad), np.ravel(w_tan)




def M_tab(pa,inc,x0,y0,theta_b,rings, delta,k, shape, vel_map, pixel_scale=1,vmode = "radial"):

	[ny,nx] = shape


	pa, inc , phi_b = pa*np.pi/180, inc*np.pi/180 , theta_b*np.pi/180 
	X = np.arange(0, nx, 1)
	Y = np.arange(0, ny, 1)
	xy_mesh = np.meshgrid(X,Y)


	e_vel = np.ones((ny,nx))*5
	vel_val = vel_map

	if vmode == "circular":


		weigths_k,weigths_j,mask = weigths_w(xy_mesh,shape,pa,inc,x0,y0,rings,delta,k,pixel_scale =pixel_scale)
		weigths_k,weigths_j = np.asarray(weigths_k),np.asarray(weigths_j)


		w_rot = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0,mask,vmode)


		sigma_v = e_vel[mask]
		x11,x12 = w_rot**2/sigma_v**2,0/sigma_v**2
		x21,x22 = 0/sigma_v**2,0/sigma_v**2


		D = (vel_val[mask])
		y1 = (w_rot/sigma_v**2)*D
		y2 = (0/sigma_v**2)*D

		A = np.asarray([[np.nansum(x11),np.nansum(x12)],[np.nansum(x21),np.nansum(x22)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2)])


		vrot,vrad = np.nansum(y1)/np.nansum(x11), 0

		
		return vrot,0,0
		#"""



	if vmode == "radial":

		weigths_k,weigths_j,mask = weigths_w(xy_mesh,shape,pa,inc,x0,y0,rings,delta,k,pixel_scale =pixel_scale)
		weigths_k,weigths_j = np.asarray(weigths_k),np.asarray(weigths_j)


		w_rot, w_rad = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0,mask,vmode)



		sigma_v = e_vel[mask]
		x11,x12 = w_rot**2/sigma_v**2,w_rot*w_rad/sigma_v**2
		x21,x22 = w_rot*w_rad/sigma_v**2,w_rad**2/sigma_v**2


		D = (vel_val[mask])
		y1 = (w_rot/sigma_v**2)*D
		y2 = (w_rad/sigma_v**2)*D

			

		A = np.asarray([[np.nansum(x11),np.nansum(x12)],[np.nansum(x21),np.nansum(x22)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2)])


		x = np.linalg.solve(A, B)
		vrot,vrad = abs(x[0]),x[1]


		if np.isfinite(vrot) == False: vrot = 0
		if np.isfinite(vrad) == False: vrad = 0

		return vrot,vrad,0

	if vmode == "bisymmetric":

		weigths_k,weigths_j,mask = weigths_w(xy_mesh,shape,pa,inc,x0,y0,rings,delta,k,pixel_scale =pixel_scale)
		weigths_k,weigths_j = np.asarray(weigths_k),np.asarray(weigths_j)


		w_rot, w_rad, w_tan = trigonometric_weights(xy_mesh,pa,inc,x0,y0,0,mask,vmode)



		sigma_v = e_vel[mask]
		x11,x12,x13 = w_rot**2/sigma_v**2,w_rot*w_rad/sigma_v**2,w_tan*w_rot/sigma_v**2
		x21,x22,x23 = w_rot*w_rad/sigma_v**2,w_rad**2/sigma_v**2,w_rad*w_tan/sigma_v**2
		x31,x32,x33 = w_rot*w_tan/sigma_v**2,w_rad*w_tan/sigma_v**2,w_tan**2/sigma_v**2


		D = (vel_val[mask])
		y1 = (w_rot/sigma_v**2)*D
		y2 = (w_rad/sigma_v**2)*D
		y3 = (w_tan/sigma_v**2)*D


		A = np.asarray([[np.nansum(x11),np.nansum(x12),np.nansum(x13)],[np.nansum(x21),np.nansum(x22),np.nansum(x23)],[np.nansum(x31),np.nansum(x32),np.nansum(x33)]])
		B= np.asarray([np.nansum(y1),np.nansum(y2),np.nansum(y3)])



		try:
			x = np.linalg.solve(A, B)
			vrot,vrad,vtan = abs(x[0]),x[1],x[2]
			if np.isfinite(vrot) == False: vrot = 50
			if np.isfinite(vrad) == False: vrad = 0
			if np.isfinite(vtan) == False: vtan = 0

		except(TypeError):
			w_sys,w_rot,w_rad,w_tan,vrot,vrad,vsys,vtan =  0,0,0,0,0,0,0,0

		return vrot,vrad,vtan



