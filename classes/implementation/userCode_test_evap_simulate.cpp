double naval = sharedParamNum(na_val);
double i_radex = inputExt(radex);
double i_glorad_max = inputExt(glorad_max);
double i_longrad = inputExt(rad_long);

// account for sub-daily resolution in calculation of energy budget (i.e. cloudiness correction factor); only needed for long-wave radiation calculation
if ( (abs(i_longrad - naval) < 0.01) && (abs(inputExt(rad_net) - naval) < 0.01) ) {
	// calculate short-wave radiation under clear sky if not given
	if ( abs(i_glorad_max - naval) < 0.01 ) {
		// calculate radex if not given
		if ( abs(i_radex - naval) < 0.01 )
			i_radex = rad_extraterr_hourly(inputExt(doy), paramNum(lat), inputExt(hour), inputExt(utc_add), paramNum(lon));
		
		i_glorad_max = calc_glorad_max(sharedParamNum(choice_gloradmax), i_radex, sharedParamNum(radex_a), sharedParamNum(radex_b), paramNum(elev));
	}
	// if i_glorad_max is very small assume nighttime and do not update state of longrad (i.e. assume persistance of cloudiness correction factor over nighttime)
	// take 5 W/m2 as "small" (although this is somewhat arbitrary) to avoid impact of atmospheric disturbances during low sun angle
	if ( i_glorad_max > 5. ) {
		i_longrad = net_longrad(inputExt(temper),inputExt(rhum),inputExt(glorad),i_glorad_max,
																					 sharedParamNum(emis_a),sharedParamNum(emis_b),
																					 sharedParamNum(fcorr_a),sharedParamNum(fcorr_b));
	} else
		i_longrad = stateScal(s_longrad);
}

set_stateScal(s_longrad) = i_longrad;



// potential evap in m/s
double et_p = et_pot (
// Common meteorological variables
	inputExt(temper),									// Air temperature (°C)
	inputExt(temp_min),								// Minimum temperature within time step (°C)
	inputExt(temp_max),								// Maximum temperature within time step (°C)
  inputExt(glorad),   						// Downward short-wave radiation (W/m2)
	inputExt(rhum),										// relative humidity (%)
	inputExt(wind),										// wind speed (m/s)
	inputExt(apress),       						// Air pressure (hPa)
	inputExt(sundur),									// Sunshine duration of current day (h)
	inputExt(cloud),										// Cloudiness (%)
// Common site-secific parameters
	paramNum(lat),											// Latitude (decimal degree)
	paramNum(lon),											// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
	paramNum(elev),										// Elevation above sea level (m)
// Specific meteorological quantities (commonly calculated internally)
	i_radex,														// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	i_glorad_max, 							// Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	inputExt(rad_net),										// net incoming ( (1-alb) * short-wave + long-wave) radiation (Wm-2)
	inputExt(rad_net_soil),									// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	stateScal(s_longrad),									// Net incoming long-wave radiation (Wm-2)
	inputExt(soilheat),								// Soil heat flux (Wm-2)
	inputExt(totalheat),								// heat conduction into soil AND plants due to physical and biochemical energy storage (Wm-2)
// Meteorological parameters
	sharedParamNum(h_tempMeas),							// height of temperature measurement above ground (m)
	sharedParamNum(h_humMeas),								// height of humidity measurement above ground (psychrometer) (m)
	sharedParamNum(h_windMeas),							// height of windspeed measurement above ground (m)
	sharedParamNum(emis_a),									// Coefficient a for calculating net emissivity (-)
	sharedParamNum(emis_b),									// Coefficient b for calculating net emissivity (-)
	sharedParamNum(fcorr_a),									// Coefficient a for calculating cloudiness correction factor (-)
	sharedParamNum(fcorr_b),									// Coefficient b for calculating cloudiness correction factor (-)
	sharedParamNum(radex_a),									// Angstrom coefficient (share of radex on glorad under clouds) (-)
	sharedParamNum(radex_b),									// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	sharedParamNum(f_day),										// soil heat flux calculation: Fraction of net_rad over daytime (in case of sub-daily application) (-)
	sharedParamNum(f_night),									// soil heat flux calculation: Fraction of net_rad over nighttime (in case of sub-daily application) (-)
// Vegetation and land-cover parameters and variables
	paramNum(crop_makk),     					// Crop-factor after Makkink (-)
	paramNum(crop_faoref),     				// Crop-factor for FAO reference approach (-)
	inputExt(cano_height),							// canopy height (m)
	inputExt(lai),											// leaf area index (m2/m2)
	inputExt(alb),											// albedo (-)
	sharedParamNum(ext),											// Canopy extinction coefficient, Beer's law (-) -  in original WASA model code set to 0.5
	paramNum(res_leaf_min),						// Plant-specific minimum (i.e. no stress occurs) stomatal resistance of a single leaf (sm-1)
	paramNum(soil_dens),								// bulk density of soil of topmost soil horizon (kg/m3)
	paramNum(glo_half),								// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA model a value of 100
	sharedParamNum(res_b),										// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
	sharedParamNum(drag_coef),								// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	sharedParamNum(rough_bare),							// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	sharedParamNum(eddy_decay),							// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	sharedParamNum(rss_a),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to 26
	sharedParamNum(rss_b),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to -1
// Computational parameters
	inputExt(doy),											// Current day of year
	inputExt(hour),												// hour of day in local time (including daylight daving time), range of [0..23]
	inputExt(utc_add),										// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	naval,															// NA value
	delta_t,						// time step length (s)
// Choice flags
	sharedParamNum(choice_et),											// flag which method to use
	sharedParamNum(choice_rcs),											// Flag: canopy stomatal resistance by up-scaling of single leaf resistance, 1: Shuttlewort & Wallace (1985) eq. 19, 2: Saugier and Katerji (1991) eq. 4 used in WASA
	sharedParamNum(choice_roughLen), 								// Flag: Roughness length for momentum transfer, 1: SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5, 2: Shuttleworth & Gurney (1990) eq. 43
	sharedParamNum(choice_plantDispl),							// Flag: Displacement height for a plant, 1: SWAT manual (2011) eq. 2:2.2.7, 2: Shuttleworth & Gurney (1990) eq. 42
	sharedParamNum(choice_gloradmax)								// Flag: calculation of maximum incoming short-wave radiation (clear sky), 1: Angstroem, 2: Allen (2005), ASCE standard etp, eq. 19 (based on elevation)
);

// output as sum in mm
set_output(etp) = et_p * delta_t * 1000.;





// potential evap in m/s
double et_a = et_act (
// Parameters and variables needed specifically for et_act()
	inputExt(wc_vol_top),							// Actual volumetric soil water content at topmost horizon (m3/m3)
	inputExt(wc_vol_root),							// Actual volumetric soil water content of the root zone (m3/m3)
	paramNum(wc_sat),									// Volumetric water content at saturation of the root zone (m3/m3)
	paramNum(wc_pwp),									// Volumetric water content of the root zone at permanent wilting point (m3/m3)
	paramNum(wc_res),									// Residual volumetric water content of the root zone (m3/m3)
	paramNum(wc_etmax),								// Parameter giving the volumetric water content where et_act equals et_pot, typically wc_etmax / wc_fk = [0.5..0.8] (m3/m3)
	paramNum(bubble),									// Bubbling pressure of the root zone (cm) or (hPa)
	paramNum(pores_ind),								// Pore-size-index of the root zone (-)
	paramNum(wstressmin),							// Threshold for water stress effect on resistance (begin of stomata closure) (cm) OR (hPa)
	paramNum(wstressmax),							// Threshold for water stress effect on resistance (total stomata closure, wilting point) (cm) OR (hPa)
	paramNum(par_stressHum),						// Parameter to calculate water vapour deficit stomatal conductance stress factor (hPa-1) - in WASA a value of 0.03
// Common meteorological variables
	inputExt(temper),									// Air temperature (°C)
	inputExt(temp_min),								// Minimum temperature within time step (°C)
	inputExt(temp_max),								// Maximum temperature within time step (°C)
  inputExt(glorad),   						// Downward short-wave radiation (W/m2)
	inputExt(rhum),										// relative humidity (%)
	inputExt(wind),										// wind speed (m/s)
	inputExt(apress),       						// Air pressure (hPa)
	inputExt(sundur),									// Sunshine duration of current day (h)
	inputExt(cloud),										// Cloudiness (%)
// Common site-secific parameters
	paramNum(lat),											// Latitude (decimal degree)
	paramNum(lon),											// Longitude of the location of interest (decimal degrees west of Greenwich, e.g. Greenwich: 0°, Berlin: 345°, New York: 75°)
	paramNum(elev),										// Elevation above sea level (m)
// Specific meteorological quantities (commonly calculated internally)
	i_radex,														// Incoming extraterrestrial radiation (i.e. at top of atmosphere) (W/m2)
	i_glorad_max, 							// Downward short-wave radiation under clear (cloudless) sky (Wm-2)
	inputExt(rad_net),										// net incoming (short-wave + long-wave) radiation (Wm-2)
	inputExt(rad_net_soil),									// net incoming (short-wave + long-wave) radiation hitting the soil surface (Wm-2)
	stateScal(s_longrad),									// Net incoming long-wave radiation (Wm-2)
	inputExt(soilheat),								// Soil heat flux (Wm-2)
	inputExt(totalheat),								// heat conduction into soil AND plants due to physical and biochemical energy storage (Wm-2)
// Meteorological parameters
	sharedParamNum(h_tempMeas),							// height of temperature measurement above ground (m)
	sharedParamNum(h_humMeas),								// height of humidity measurement above ground (psychrometer) (m)
	sharedParamNum(h_windMeas),							// height of windspeed measurement above ground (m)
	sharedParamNum(emis_a),									// Coefficient a for calculating net emissivity (-)
	sharedParamNum(emis_b),									// Coefficient b for calculating net emissivity (-)
	sharedParamNum(fcorr_a),									// Coefficient a for calculating cloudiness correction factor (-)
	sharedParamNum(fcorr_b),									// Coefficient b for calculating cloudiness correction factor (-)
	sharedParamNum(radex_a),									// Angstrom coefficient (share of radex on glorad under clouds) (-)
	sharedParamNum(radex_b),									// Angstrom coefficient (radex_a+radex_b = share of radex on glorad under clear sky) (-)
	sharedParamNum(f_day),										// soil heat flux calculation: Fraction of net_rad over daytime (in case of sub-daily application) (-)
	sharedParamNum(f_night),									// soil heat flux calculation: Fraction of net_rad over nighttime (in case of sub-daily application) (-)
// Vegetation and land-cover parameters and variables
	paramNum(crop_makk),     					// Crop-factor after Makkink (-)
	paramNum(crop_faoref),     				// Crop-factor for FAO reference approach (-)
	inputExt(cano_height),							// canopy height (m)
	inputExt(lai),											// leaf area index (m2/m2)
	inputExt(alb),											// albedo (-)
	sharedParamNum(ext),											// Canopy extinction coefficient, Beer's law (-) -  in original WASA model code set to 0.5
	paramNum(res_leaf_min),						// Plant-specific minimum (i.e. no stress occurs) stomatal resistance of a single leaf (sm-1)
	paramNum(soil_dens),								// bulk density of soil of topmost soil horizon (kg/m3)
	paramNum(glo_half),								// Solar radiation at which stomatal conductance is half of its maximum (W/m2) - in WASA model a value of 100
	sharedParamNum(res_b),										// Mean boundary layer resistance (sm-1) - in SW and WASA a value of 25
	sharedParamNum(drag_coef),								// Effective value of mean drag coef. of vegetative elements (-) - in SG and WASA a value of 0.07
	sharedParamNum(rough_bare),							// Roughness length of bare substrate (m) - in SW, SG and WASA a value of 0.01
	sharedParamNum(eddy_decay),							// Eddy diffusivity decay constant (-) - by SW and SG set to 2.5 for agricultural crops
	sharedParamNum(rss_a),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to 26
	sharedParamNum(rss_b),										// Empirical constant for calculation of soil surface resistance (-) - in original code set to -1
// Computational parameters
	inputExt(doy),											// Current day of year
	inputExt(hour),												// hour of day in local time (including daylight daving time), range of [0..23]
	inputExt(utc_add),										// Deviation of local time zone from UTC; may vary over the year due to daylight saving time; range of [-12..14] (hours)
	naval,														// NA value
	delta_t,						// time step length (s)
// Choice flags
	sharedParamNum(choice_et),											// flag which method to use
	sharedParamNum(choice_rcs),											// Flag: canopy stomatal resistance by up-scaling of single leaf resistance, 1: Shuttlewort & Wallace (1985) eq. 19, 2: Saugier and Katerji (1991) eq. 4 used in WASA
	sharedParamNum(choice_roughLen), 								// Flag: Roughness length for momentum transfer, 1: SWAT manual (2011) eqs. 2:2.2.4, 2:2.2.5, 2: Shuttleworth & Gurney (1990) eq. 43
	sharedParamNum(choice_plantDispl),							// Flag: Displacement height for a plant, 1: SWAT manual (2011) eq. 2:2.2.7, 2: Shuttleworth & Gurney (1990) eq. 42
	sharedParamNum(choice_gloradmax)								// Flag: calculation of maximum incoming short-wave radiation (clear sky), 1: Angstroem, 2: Allen (2005), ASCE standard etp, eq. 19 (based on elevation)
);

// output as sum in mm
set_output(eta) = et_a * delta_t * 1000.;
