/////////////////////////////////////////////////////////////////////////////////
// Authors: Tobias Pilz, Julius Eberhard (comments)
// Last Edit: Unknown
// Project: ECHSE Evapotranspiration
// Engine: soilheat_portugal
// Method: simulate
// Aim: Calculating Soil Heat Flux
/////////////////////////////////////////////////////////////////////////////////

// calculate radex
double i_radex = rad_extraterr_hourly(inputExt(doy),
                                      paramNum(lat),
                                      inputExt(hour),
                                      inputExt(utc_add),
                                      paramNum(lon));

// calculate glorad_max
double i_glorad_max = calc_glorad_max(sharedParamNum(choice_gloradmax),
                                      i_radex,
                                      sharedParamNum(radex_a),
                                      sharedParamNum(radex_b),
                                      paramNum(elev));

// determine value of i_daynight
double i_daynight = -9999.;
if (i_glorad_max > 1e-6)
	i_daynight = 1;
else
	i_daynight = 0;

// net radiation at soil surface
double i_rad_net = -9999.;
if (sharedParamNum(choice_et) == 13)
	// canopy extinction according to Beer's law:
	// net radiation at soil surface (W.m-2)
  i_rad_net = inputExt(rad_net) * exp(-1. * sharedParamNum(ext) * inputExt(lai));
else
  i_rad_net = inputExt(rad_net);

// calculate soil heat flux
double i_soilheat = soil_heatflux(i_rad_net,
                                  sharedParamNum(f_day),
                                  sharedParamNum(f_night),
                                  i_daynight,
                                  delta_t);

set_output(soilheat_out) = i_soilheat;
