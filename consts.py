# https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
# M_sat = 5.6834e26
# M_jup = 1.89812e27
M_sun = 1.98841e30

# masses relative to the sun
M_sun = 1.00000597682 # M_Sun normalized to 1
M_jup = 0.000954786104043 # M_Jup / M_Sun
M_sat = 0.000285583733151 # M_Sat / M_Sun
M_earth = 3.003e-6 # M_earth / M_Sun

# 1 m = 6.6845e-12 au
meter_to_au = 6.6845871222684e-12
# 1au = 1.495e11 m
au_to_meter = 149597870700

# 1 second = 1 / 86400 day
second_to_day = 1 / 86400
# 1 day = 24 * 60 * 60 s
day_to_second = 86400

# [G]: m^3 / kg s^2 => using a.u and days [G]: a.u / kg days^2
# G = (6.67430e-11 * (meter_to_au ** 3) * (day_to_second ** 2))

# [G] = au^3 / ( (kg / M_Sun) * day^2)
G = 2.95912208286e-4

# initial conditions
# positions are in a.u (astronomical units)
# speeds are in a.u / day => speed * mass = impulsion
# masses are in kg
# http://vo.imcce.fr/webservices/miriade/?forms

sun_position0 = [0., 0., 0.]
sun_impulsion0 = [0., 0., 0.]

jupiter_position0 = [3.4707903364632, -3.3666298150704, -1.5275164476857]
jupiter_impulsion0 = [0.0054161298800 * M_jup, 0.0051281500076 * M_jup, 0.0020662323714 * M_jup]

saturn_position0 = [5.8139930169916, -7.3998049325634, -3.3068297277913]
saturn_impulsion0 = [0.0042287765130 * M_sat, 0.0030656447687 * M_sat, 0.0010842701680 * M_sat]

earth_position0 = [ -0.6328481424305, -0.7224546684653, -0.3131787271070]
earth_impulsion0 = [0.0131352571260 * M_earth, -0.0099512546159 * M_earth, -0.0043145281044 * M_earth]

# time in days
t0 = 0
tN = 5000 * 365.25 # integration over 5000 years
dt = 30 # time_step : 30 days

# refresh times for animation
DATA_PLOT_REFRESH = 30 # 30 days
DATA_SUB_INTERVAL_LENGTH = 12 * DATA_PLOT_REFRESH # simulation subinterbal is ~ 1 year length
