from readNC import *

varNames = ["prec","evap", "runoff", "baseflow", "wdew", "sm1", "sm2", "sm3", "evap_canop", "evap_veg", "evap_bare", "net_short", "net_long", "r_net", "surf_temp", "swq", "snow_depth", "snow_canop"]
fileName = 'PGF/VIC/output/output_grid_1981010100.ctl'
data = readGrads(fileName, "prec", str(1), lon=[-179.5, 179.5])
createNetCDF("output_19810101.nc", varNames, ["mm","mm","mm","mm","mm","mm","mm","mm","mm","mm","mm","W/m^2","W/m^2","W/m^2", "C", "mm","cm","mm"], latitudes=data.latitudes, longitudes=data.longitudes, loop=True)

for i in range(365):
    print i
    for var in varNames:
       data = readGrads(fileName, var, str(i+1), lon=[-179.5, 179.5])
       data2NetCDF("output_19810101.nc", var, data, data.time, posCnt = i)

