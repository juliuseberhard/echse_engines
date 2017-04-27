/////////////////////////////////////////////////////////////////////////////////
// Author: Julius Eberhard
// Last Edit: 2016-10-12
// Project: ECHSE Evapotranspiration
// Engine: glorad_portugal
// Method: simulate
// Aim: Calculating Global Radiation
/////////////////////////////////////////////////////////////////////////////////

// calculate extraterrestrial radiation daily (not hourly!)
double i_radex = rad_extraterr_daily(inputExt(doy),
                                     paramNum(lat));

// calculate global radiation
double i_glorad = calc_glorad(i_radex,
                              inputExt(sundur),
                              inputExt(cloud),
                              paramNum(lat),
                              inputExt(doy),
                              sharedParamNum(radex_a),
                              sharedParamNum(radex_b));

// output
set_output(glorad_out) = i_glorad;
