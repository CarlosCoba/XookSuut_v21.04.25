import numpy as np
from astropy.io import fits



def save_vlos_model(galaxy,vmode,vel_map,evel_map,model2D,PA,INC,XC,YC,VSYS,save = 1):

	if save == 1:

		# The best 2D model
		hdu = fits.ImageHDU()
		hdu.data = model2D

	
		hdu.header['NAME0'] = '2D LoS velocity model'
    
		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s.vlos_model.fits"%(galaxy,vmode),overwrite=True)


		
		# Now the residual map

		hdu = fits.ImageHDU()
		hdu.data = vel_map- model2D

	
		hdu.header['NAME0'] = 'residual map, data-model'
    
		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s.residual.fits"%(galaxy,vmode),overwrite=True)


		
		# Now the chisquare map

		sigma = evel_map
		chisq = ( vel_map- model2D )/sigma
		chisq_2 = chisq**2
		hdu = fits.ImageHDU()
		hdu.data = chisq_2

	
		hdu.header['NAME0'] = 'Chisquare map, (data-model)**2/sigma**2'
    
		hdu.header['PA'] = PA
		hdu.header['INC'] = INC
		hdu.header['VSYS'] = VSYS
		hdu.header['XC'] = XC
		hdu.header['YC'] = YC
		
		hdu.writeto("./fits/%s.%s.chisq.fits"%(galaxy,vmode),overwrite=True)



	else: pass

