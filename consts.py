from astropy.constants import M_sun, M_jup

# https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
M_sat = 5.6834e26

# initial conditions 
# positions are in a.u (astronomical units)
# speeds are in a.u / day => speed * mass = impulsion
# masses are in kg
# http://vo.imcce.fr/webservices/miriade/?forms

sun_position0 = [0, 0, 0]
sun_impulsion0 = [0, 0, 0]

jupiter_position0 = [3.4707903364632, -3.3666298150704, -1.5275164476857]
jupiter_impulsion0 = [0.0054161298800 * M_jup, 0.0051281500076 * M_jup, 0.0020662323714 * M_jup]

saturn_position0 = [5.8139930169916, -7.3998049325634, -3.3068297277913]
saturn_impulsion0 = [0.0042287765130 * M_sat, 0.0030656447687 * M_sat, 0.0010842701680 * M_sat]

t0 = 0
tN = 5000
dt = 1