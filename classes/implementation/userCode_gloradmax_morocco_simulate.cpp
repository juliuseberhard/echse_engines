/////////////////////////////////////////////////////////////////////////////////
// Author: Julius Eberhard
// Last Edit: 2016-11-28
// Project: ECHSE Evapotranspiration
// Engine: gloradmax_morocco
// Method: simulate
// Aim: Calculating Clear Sky Radiation
/////////////////////////////////////////////////////////////////////////////////

double glorad_max = -9999.;

// calculate extraterrestrial radiation hourly
double i_radex = rad_extraterr_hourly(inputExt(doy),
                                      paramNum(lat),
                                      inputExt(hour),
                                      inputExt(utc_add),
                                      paramNum(lon));

// calculate clear sky global radiation
double i_glorad_max = calc_glorad_max(sharedParamNum(choice_gloradmax),
                                      i_radex,
                                      sharedParamNum(radex_a),
                                      sharedParamNum(radex_b),
                                      paramNum(elev));

// update glorad_max, if it is small (< 20 W.m-2) during nighttime
if ( i_glorad_max < 20. ) {
  glorad_max = stateScal(s_glorad_max);
} else {
	glorad_max = i_glorad_max;
}

set_output(gloradmax_out) = glorad_max;
