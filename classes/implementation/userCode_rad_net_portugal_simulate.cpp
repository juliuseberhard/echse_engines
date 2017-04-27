/////////////////////////////////////////////////////////////////////////////////
// Author: Julius Eberhard
// Last Edit: 2016-11-28
// Project: ECHSE Evapotranspiration
// Engine: rad_net_portugal
// Method: simulate
// Aim: Calculating Net Radiation
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

set_stateScal(s_glorad_max);

// calculate net radiation
double i_rad_net = net_rad(inputExt(glorad),
                           glorad_max,
                           inputExt(temper),
                           inputExt(rhum),
                           inputExt(alb),
                           sharedParamNum(emis_a),
                           sharedParamNum(emis_b),
                           sharedParamNum(fcorr_a),
                           sharedParamNum(fcorr_b));

set_output(rad_net_out) = i_rad_net;
