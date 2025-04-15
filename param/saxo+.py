#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Isaac Dinis
"""



import shesha.config as conf
import numpy as np

simul_name = "saxo+"  

fs = 3000 # [Hz]


# loop
p_loop = conf.ParamLoop()
p_loop.set_niter(5000)                  
p_loop.set_ittime(1./freq_second_stage)
p_loop.set_devices([0])      

# geom
p_geom = conf.ParamGeom()
p_geom.set_pupdiam(400)         # /!\ forcing even when the number of pixels of the pyr is reduced
p_geom.set_zenithangle(0.)      # /!\ keep it 0 until we know what it does in COMPASS

# tel
p_tel = conf.ParamTel()
p_tel.set_diam(8.0)            # /!\  VLT diameter
p_tel.set_cobs(0.14)           # /!\  central obstruction
p_tel.set_type_ap("VLT")       # /!\  VLT pupil
p_tel.set_t_spiders(0.00625)   # /!\  spider width = 5 cm

# atmos
p_atmos = conf.ParamAtmos()   # /!\ Dummy atmosphere for setup of WFS ...
p_atmos.set_r0(0.14)          # /!\ Fried parameter will be overwritten by the SAXO.PY value
p_atmos.set_nscreens(1)       # /!\ Number of layers
p_atmos.set_frac([1.0])       # /!\ Fraction of atmosphere (100% = 1)
p_atmos.set_alt([0.0])        # /!\ Altitude(s) in meters
p_atmos.set_windspeed([10])   # /!\ Wind speed of layer(s) in m/s
p_atmos.set_winddir([45])     # /!\ Wind direction in degrees
p_atmos.set_L0([25])          # /!\ Outer scale in meters

# target
p_target = conf.ParamTarget()
p_targets = [p_target]
p_target.set_xpos(0.)         # /!\ On axis
p_target.set_ypos(0.)         # /!\ On axis
p_target.set_Lambda(1.65)     # [microns] H band.
p_target.set_mag(6.)          # /!\

# wfs
p_wfs0 = conf.ParamWfs(roket=True)
p_wfss = [p_wfs0]

p_wfs0.set_type("pyrhr")        # /!\ pyramid
p_wfs0.set_nxsub(50)            #     number of pixels
p_wfs0.set_fracsub(0.0001)      #     threshold on illumination fraction for valid pixel
p_wfs0.set_Lambda(1.043)        # [microns] Y bandwidth (width = 140nm)

p_wfs0.set_gsmag(6.)            
p_wfs0.set_zerop(1.e11)         
p_wfs0.set_optthroughput(1)     

# Standard deviation of additive noise
# If <0, no noise, if =0 photon noise enabled, if >0 add additive noise
# When a (positive) readout noise is set here, the contributions of the sky 
# background and the dark current will be added, taking duty cycle and excess 
# noise factor of the camera. This will be added after the actual frequency of 
# the camera is known.
p_wfs0.set_noise(simu.parameter("stage2_readout_noise"))
# p_wfs0.set_noise(-1)         # EMCCD with < 0.1e- RON
p_wfs0.set_xpos(0.)             # /!\ On axis
p_wfs0.set_ypos(0.)             # /!\ On axis
rMod = simu.parameter("stage2_modulationRadius") # Modulation radius, in lam/D units
p_wfs0.set_pyr_ampl(rMod)
if (rMod==0):                    # generalize parameter file for whatever modulation radius value.
    nbPtMod = 1
else:
    nbPtMod = int(np.ceil(int(rMod * 2 * 3.141592653589793) / 4.) * 4)
p_wfs0.set_pyr_npts(nbPtMod)
p_wfs0.set_pyr_pup_sep(32)      # lateral shift of pupils centers (50/2+14/2 pixels for C-red one).
p_wfs0.set_fstop("round")
p_wfs0.set_fssize(3.2)          # Size of posterior field stop (arcseconds) accord. to optical design
#p_wfs0.pyr_loc("after")        # modulator position w.r.t the field stop: ["before","after"]
                                # Problem: COMPASS code always simulated "before".
p_wfs0.set_pyr_compute_focalplane(False) # in case one wants to check the PSF on the pyramid
p_wfs0.set_atmos_seen(1)        # /!\

# dm
p_dm0 = conf.ParamDm()
p_dms = [p_dm0]

p_dm0.set_type("pzt")           # /!\
p_dm0.set_thresh(-200)          # /!\ to get all Boston actuators
p_dm0.set_alt(0.)               # /!\
p_dm0.set_unitpervolt(1.)       # /!\
p_dm0.set_push4imat(1.0e-3)
p_dm0.set_file_influ_fits(os.path.join(_my_dm_dir, boston_dm_file)) # /!\ choice made at the begin. of this file
boston_dms = ['Boston24x24_flat.fits',
              'Boston28x28_flat.fits',
              'Boston34x34_flat.fits',
              'Boston34x34_measured_IF.fits']
if (boston_dm_file in boston_dms) and recommended_dm_pup:
    # recommanded pupil diam on DMs
    boston_rec_diam = [(24-3)*450e-6, (28-3)*450e-6,
                       (34-3)*400e-6, (34-3)*400e-6]
    p_dm0.set_diam_dm(boston_rec_diam[boston_dms.index(boston_dm_file)])
else:
    print('No recommended pupil size used on Boston DM')

# centroiders
p_centroider0 = conf.ParamCentroider()
p_centroiders = [p_centroider0]

p_centroider0.set_nwfs(0)           # /!\
p_centroider0.set_type("maskedpix")

# controllers
p_controller0 = conf.ParamController()
p_controllers = [p_controller0]

p_controller0.set_type("generic")   # /?\ generic => must do manual imat.
p_controller0.set_nwfs([0])         # /!\
p_controller0.set_ndm([0])          # /!\
p_controller0.set_delay(saxoplus_compute_delay(freq_second_stage)) # /?\ needs to be refined [?]
p_controller0.set_gain(simu.parameter("stage2_gain"))

# coronagraph
p_corono0 = conf.ParamCoronagraph()
p_coronos = [p_corono0]

p_corono0.set_type(simu.parameter("coro_type"))
p_corono0.set_dim_image(200)    # size of the science image in pixel
p_corono0.set_wavelength_0(simu.parameter("lambda_science"))
                                # coronagraph wavelength in microns
                                # 1.667 μm = central wavelength of H3 IRDIS filter
                                # this value is just an example, not a fixed parameter

# this calculation is done to match IRDIS pixel scale, according to wavelength_0
IRDIS_pixel_scale = 12.25 / 1000 / 3600  # [rad] IRDIS pixel scale
VLT_diameter = 8                         # [m]
wavelength_0 = p_corono0.get_wavelength_0() * 1e-6
image_sampling = np.rad2deg(wavelength_0 / VLT_diameter) / IRDIS_pixel_scale
p_corono0.set_image_sampling(image_sampling)  # size of lambda/D in pixel

# for a polychromatic coronagraph (optional), set these two parameters:
# p_corono0.set_nb_wav(3)           # number of simulated wavelength
# p_corono0.set_delta_wav(0.054)    # spectral bandwidth in micron
                                    # 0.054 μm = bandwidth of H3 IRDIS filter
                                    # this value is just an example, not a fixed parameter
