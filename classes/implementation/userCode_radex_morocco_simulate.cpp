/////////////////////////////////////////////////////////////////////////////////
// Author: Julius Eberhard
// Last Edit: 2016-10-11
// Project: ECHSE Evapotranspiration
// Engine: radex_morocco
// Method: simulate
// Aim: Calculating Extraterrestrial Radiation
/////////////////////////////////////////////////////////////////////////////////

double i_radex = rad_extraterr_hourly(inputExt(doy),
                                      paramNum(lat),
                                      inputExt(hour),
                                      inputExt(utc_add),
                                      paramNum(lon));

set_output(radex_out) = i_radex;
