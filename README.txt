
Installation

(1) Make a directory to put the executables in:
mkdir ~/XS-code

(2) copy all files in this new folder

(3) Make the file executable:
chmod +x XookSuut.py


(4) Edit  .bashrc:
gedit ~/.bashrc

Add the folowing line:
export PATH="$PATH:$HOME/XS-code"

Refresh the bash file:
source ~/.bashrc

Try it. Go to any folder and type  XookSuut

you must get the following :


USE: XookSuut.py object vel_map [evel_map] SN [VSYS] PA INC X0 Y0 [N_it=5] [pixel_scale=1] vary_PA[0,1] vary_INC[0,1] vary_XC[0,1] vary_YC[0,1] vary_VSYS[0,1] vary_THETA_bar[0,1] ring_space [Delta] [Rstart=2.5] [Rfinal=40] [frac_pixel=2/3.] [R_bar_min,R_bar_max] [model=cicular,radial,bisymmetric] [save_plots=1] [errors=0] [survey] [config_file] [e_ISM=5]")



"""
#################################################################################
# 				The XookSuut-code.				#
#				version 21.04.25				#
# 				C. Lopez-Coba					#
#################################################################################


INPUT:
-----------------------------------------------------------------
parameter|	default value	|	Description 
-----------------------------------------------------------------
object =	- 		Object name.
vel_map = 	-		2D velocity map in km/s.
evel_map = 	-		2D error map in km/s. If no error map is passed then set this to "".
SN = 		-		S/N cut applied on the velocity map. 
				Better results are obtained when SN has low values, 15-20 km/s.
				if evel_map is set to 0, then the whole velocity map will be used for estimating the rotation curve.
VSYS	=	""		Guess value for Vsys. If VSYS is set to "", then it will be used the central value as guess (i.e., VSYS =  vel_map[Y0][X0]).
PA	= 	-		Guess value for the position angle of the major axis.		
INC	=	-		Guess value for the inclination.
X0,Y0	=	-		Guess values for the kinematic center. These are in pixel coordinates. 
N_it 	=	5		Number of iteratios to find the best parameters. 
pixel_scale	1		Size of the pixel in arcsecs.
vary_PA =	1		Boolean.
vary_INC =	1		Boolean.
vary_XC =	1		Boolean.
vary_YC =	1		Boolean.
vary_VSYS =	1		Boolean.
vary_PHI_BAR =	1		Boolean. Position angle of the bisimetric distortion on the galaxy plane.
delta =		""		Width of the rings (2*delta) in the tabulated model. If delta is set to "", then delta = 0.5*ring_space
ring_space =	2.5		Spacing between rings used to create the interpolated model. 
rstart	=	2.5		Initial radius in the interpolated model. 
rfinal	=	40		Final radius in the interpolated model.
frac_pixel = 	2/3.		Fraction of total pixels within a ring required to compute the tabulated model. If frac_pixel = 1, 
				then 100% of the pixels within a ring are required to be filled with data to compute the tabulated model.

r_bar_min,max= 	2.5,40		Maximum and minimum length of the bisymmetric perturbation. If only a value is passed then it'll be considered as r_bar_max.
model	= 	circular	The different kinematic models. You must choose between: circular motion (circular), radial flows (radial), or barlike flows (bisymmetric).
errors	=	0		Boolean. If you want to compute errors via MCMC for the derived parameters.
survey =	-		String. If the object belongs to a specific galaxy survey.
config =	""		configure file to pass initial guess and constrains for the constant params (x0, y0, Vsys, pa, inc, phi_bar). This file is optional.
				If this file is passed to XookSuut, then it will ommit the previos guess values (VSYS, PA, INC, X0, Y0) as well as the VARY_ entrance.
save_plots =			Boolean
e_ISM =		5		Error in the emission line centroid.


A configuration file has the follow entrance

#
#
# XS config file
# Example of configuration file
# param col. are the constant parameters to fit. Do not touch this column.
# val col. corresponds to the initial values for the considered parameters
# fit col. is a boolen. If set to 1 means the parameter is fitted, other wise set 0.
# min col. is the minimum value for the considered parameter. If fit is seto to 0, the min value is not taken into account. 
# max col. is the maximum value for the considered parameter. If fit is seto to 0, the max value is not taken into account. 
#
#
#
param	val	fit	min	max
pa	35	1	0	360
inc	35	1	0	90
x0	25	1	0	50
y0	25	1	0	50
Vsys	11168.266579019075	0	0	3e6
phi_b	45	1	0	180





					OUTPUT:


-Tables:

*Table containing the best fitted values for the constant parameters (x0, y0, Vsys, pa, inc, phi_bar)
ana_kin_model.csv

*Table containing different estimations of the maximum circular velocity.
vmax_out_params.bisymmetric_model.csv


* fits files
2D array of the LoV of the best kinematic model
1D array of the different kinematic models as function of the deprojected radius


* plots
plot of the best 2D kinematic models
plot of the asymptotic velocity estimated with Vt 



************************************************
************    ADDITIONAL *********************
************************************************

The following directories must be created in the working directory.
./plots
./fits
./vmax_rturn
"""


************************************************
************    PACKAGES *********************
************************************************
python3
numpy
matplotlib
astropy
lmfit
numdifftools


#
# RUNNING THE CODE:
#
#
# EXAMPLE
#

XookSuut.py test_galaxy test.fits "" 20 "" 150 50 35 33 5 1 1 1 1 1 1 1 2.5 "" 3 46 1/3. 20 bisymmetric 1 0 test_object "" 5


