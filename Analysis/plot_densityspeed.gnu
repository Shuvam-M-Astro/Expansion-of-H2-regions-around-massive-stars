#! /usr/bin/gnuplot 
# loop through every .dat file from a simulation run and read in the data for density and speed (entries 2 and 3) as functions of radius (in # pc). Then plot both on the same figure, loglog in # density/radius and linear in speed. Ultimately to be used to compile snapshots of pngs  # into an animation through time.
# use ternary operator to append string "" or "notafile" to the .dat file to either plot ionisation fraction or not.

load("settings.gnu")

# lambda value, from bondiaccretion.pdf, eq 8.22
qs(g) = 0.25*(2./(5-3.*g))**((5.-3.*g)/(2.*g-2.)) 

# sonic/critical radius, pc
r_c(M_central, a_) = G * M_central * SolarMass / (2. * a_**2. * constant_pc)

# accretion radius, from bondiaccretion.pdf, pc
r_a(M_central, a_inf_) = 2. * G * M_central * SolarMass / (a_inf_**2. * constant_pc)

# sound speed function
soundspeed(g, T, mu) = sqrt(g * T * R_constant / mu)

# low radii approximation for flow speed, from bondiaccretion.pdf
speed_approximation(r, a_inf, r_a_) = -a_inf * (r / r_a_)**-0.5

# analytic speed relation, bondiaccretion.pdf equation 8.25. Taking sound speed constant, isothermal 
speed_analytic(g_ions_, M_central, r, a_inf_) = -qs(g_ions_) * (G * M_central * SolarMass)**2 / (r**2 * a_inf_**3)

# low radii density approximation, from bondiaccretion.pdf
density_approximation(r, r_a_, g_ions_) = 0.25 * qs(g_ions_) * rho_inf * (r / r_a_)**-1.5

# low radii density approximation, from Keto (2003) --- rho_sonicpoint is taken from 3000th frame of Bondi flow, see bondi_analyse.py
density_freefall(r, r_crit, rho_sonicpoint_) = rho_sonicpoint_ * (r / r_crit)**-1.5

# hydrostatic density profile for large radii, Keto (2003)
density_hydrostatic(r, r_crit, rho_sonicpoint_) = rho_sonicpoint_ * exp(-2 * (1 - r_crit / r))

# updated hydrostatic density profile for large radii, revision of Keto (2003) - suggestion by Rolf and Dylan
density_updated_hydrostatic(r, r_crit, rho_sonicpoint_, rho_inf_) = rho_sonicpoint_ * exp(log(rho_inf_ / rho_sonicpoint_) * (1 - r_crit / r))

# hydrostatic density profile, from MHD textbook, compare to above ^
density_hydrostatic_function2(r, r_crit, rho_sonicpoint_) = rho_sonicpoint_ * exp(r * (1 - r_crit / r))

# first order large radii density profile
density_firstorder(r, r_crit, rho_sonicpoint_) = rho_sonicpoint_ * exp(2 * r_crit / r)

# evaluate sound speed at infinity, in cm/s
a_inf = soundspeed(g_inf, T_inf, mu_inf)

# evaluate ionised gas sound speed, in cm/s
a_ions = soundspeed(g_ions, T_ions, mu_ions)

# evaluate accretion and sonic radii
r_a = r_a(M_central, a_inf)
r_c = r_c(M_central, a_ions)
r_c_neutrals = r_c(M_central, a_inf)

# sound speed at sonic point, see bondiaccretion.pdf
a_sonic(a_inf_, gamma) = a_inf_ * (2. / (5. - 3. * gamma)) ** 0.5

# sonic radius, bondiaccretion.pdf
r_sonic(a_inf_, M_central, gamma) = 0.25 * (5. - (3. * gamma)) * G * M_central * SolarMass / (a_inf_**2. * constant_pc)

# density at sonic point, bondiaccretion.pdf
rho_sonic(rho_inf_, gamma) = rho_inf * (2. / (5. - 3. * gamma)) **(1. / (gamma - 1.)) 

# Accretion rate M dot, from bondiaccretion.pdf, cgs
M_dot(g_ions_, rho_inf_, M_central, a_inf_) = 4. * pi * qs(g_ions_) * rho_inf_ * (G * M_central * SolarMass)**2. / a_inf_**3.

a_sonicpoint = a_sonic(a_inf, g_ions)

r_sonicpoint = r_sonic(a_inf, M_central, g_ions)

rho_sonic = rho_sonic(rho_inf, g_ions)

m_dot = M_dot(rho_inf, a_inf, M_central, g_ions)

# omega function, which is argument of Lambert W function to be used in initialising Bondi flow
omega(r, r_sonic_) = (r_sonic_ / r)**4. * exp(3. - 4. * r_sonic_ / r)

# analytic velocity function
velocity_analytic(r, r_sonic_, soundspeed_) = -soundspeed_ * (-lambertw(-omega(r, r_sonic_)))**0.5

# density at sonic point is scaled by ambient density
rho_sonic = rho_inf * exp(3. / 2.)

# rho function, given the velocity and calculate the value of density 
rho_analytic(r, r_sonic_, velocity, sound_speed) = - rho_sonic * r_sonic_ **2. * sound_speed / (r**2 * velocity_analytic(r, r_sonic_))  

print "sonic point of neutral flow = ", r_c_neutrals, " pc"
print "omega, argument of W(x) = ", omega(1e-2, r_c_neutrals), " ...unknown units"
print "velocity at 1e-5 pc as computed using principle branch of LambertW = ", velocity_analytic(1e-2, r_c_neutrals, a_inf), " cgs"
#print "rho at r = 1.0 pc = ", rho_analytic(1.0 * 3.086e18, r_c_neutrals * 3.086e18, velocity_analytic(1.0, r_c_neutrals, a_inf), a_inf) " cgs"
print "accretion radius = ", r_a, " pc"
print "velocity from approximation, at 1e-5 pc = ", speed_approximation(1e-5, a_inf, r_a)
print "velocity from analytic function at 1e-5 pc = ", speed_analytic(g_ions, M_central, 1e-5 * 3.086e18, a_inf), " cgs"
print "neutral sound speed = ", a_inf, " cgs"
print "ion sound speed = ", a_ions, " cgs"
print "sonic radius of ionised flow = ", r_c, " pc"
print "sound crossing time = ", 1.8775 * 3.086e18 / (a_inf * 3.15e7), " yrs"
print "analytic accretion rate = ", M_dot(g_ions, rho_inf, M_central, a_inf), " cgs"

# set up plotting commands

lw = 1.2

set terminal pngcairo dashed size 900,700 linewidth lw*2.0 font "Times,20"
set key top left
set tics scale 1.7 front

set format x "10^{%L}"
#set format x2 "10^{%L}"
set format y "10^{%L}" 
set format y2 "%g"

#set x2range [9e-6 * 206265 : 206265];
#set xrange [9e-6:1]
set xrange [xmin:xmax]

#set xtics nomirror
set ytics nomirror
set y2tics tc rgb "red"
set yrange [rho_min:rho_max]
set y2range [v_min:v_max]
#set x2tics

set mytics 10
set my2tics 4
set mxtics 10
#set mx2tics 10

if (logy == 1){ set logscale y 10 }
if (logx == 1){ set logscale x 10 }
#set logscale x2 10

set xlabel "Radius [pc]"
#set x2label "Radius [AU]"
set ylabel "Density [g cm^{-3}]"
set y2label "Velocity [km s^{-1}]" textcolo rgb "red" offset -1, 5

set rmargin 8

set key reverse Left right samplen 2
set key at graph 1.0, graph 0.7
#set key at graph 0.4, graph 0.3
set border 7 lc rgb "black"

eps = 1e-20  # gut für lw*2 mit lw = 1.2
set arrow 1 from r_c, graph eps to r_c, graph 1.0 nohead lw lw linecolor rgb "grey" dt '-' back
#set arrow 2 from r_c, graph 0.32 to r_c, graph 1.0 nohead lw lw linecolor rgb "grey" dt '-' back
set label 1 at r_c*1.09,graph 0.95 "R_{c_{ions}}" tc rgb "#999999" 

#set arrow 3 from r_a, graph eps to r_a, graph 0.18 nohead lw lw linecolor rgb "#c97446" dt '-' back
#set arrow 4 from r_a, graph 0.35 to r_a, graph 1.0 nohead lw lw linecolor rgb "#c97446" dt '-' back
#set label 2 at r_a*1.15,graph 0.89 "R_{a}" tc rgb "#c97446" 

set arrow 5 from graph 1,grap 0.0 to gr 1,gr 1 lw lw lc rgb "red" dt 1 nohead front

if (show_ions == 1) {
	
	set arrow 6 from graph 1,grap 0 to gr 1,gr 0.33 lw lw lc rgb "#239B56" dt 1 nohead front
	set label 4 at graph 1.13,graph 0.005 "Ion. Fraction, I/I_{o}" tc rgb "#239B56" rotate by 90 left
}

if (analysis == 1){

set label 5 at graph 0.54, 0.46 "Free Fall,   n_{r_{c}}(r / r_{c})^{-1.5}" font "Times, 14" tc rgb "#737373" rotate by -28 right
set label 6 at graph 0.46, 0.93 "Analytic, c_{/Symbol \245}( r / r_{a} )^{-0.5}" font "Times, 14" tc rgb "#ff471a" rotate by 28 right
#set label 6 at graph 0.62, 0.19 "Hydrostatic, n_{r_{c}}exp[2r_{c}/r]" font "Times, 14" tc rgb "#2415c1" rotate by -0.5 left
#set label 6 at graph 0.33, 0.465 "Hydrostatic," font "Times, 14" tc rgb "#2415c1" rotate by -26 left front
#set label 7 at graph 0.32, 0.43 "n_{r_{c}}exp[ln(n_{/Symbol \245}/n_{r_{c}})(1 - r_{c}/r)]" font "Times, 14" tc rgb "#2415c1" rotate by -26 left front
set label 7 at graph 0.79, 0.31 "Hydrostatic," font "Times, 14" tc rgb "#2415c1" rotate by -12 left
set label 8 at graph 0.74, 0.24 "n_{r_{c}}exp[-2(1 - r_{c}/r)]" font "Times, 14" tc rgb "#2415c1" rotate by -12 left
}

do for [i = N_min:N_max]{
  
  if (i%100==0){ print "i: ",i }
  
  file = sprintf('data.%04i.dat', i)
  outfile = sprintf('%02i.png', i)
  set output outfile

  yearnumber = sprintf('%d',i)
  #totalyearnumber = 66600 + (yearnumber - 666)
  #print "yearnumber = ", yearnumber
  
  titlename = sprintf('Time: %d yr', 100 * i)
  
  set label 9 titlename at graph 0.67, 0.73 #0.08, graph 0.33

  p file u ($1):($3 * 9.778e10/1e5) axes x1y2 with lines lw lw lc rgb "#ff3300" title 'Velocity',\
    file u ($1):($2 * 5.210e-21) axes x1y1 with lines lw lw lc rgb "black" title 'Density',\
    (show_ions?"":"notafile").file u ($1):($11 * ion_scale - ion_offset) axes x1y2 with lines lw lw lc rgb "#239B56" title 'Ion. Fraction',\
    (analysis?1:1/0) * speed_approximation(x, a_inf, r_a)/1e5 axes x1y2 with lines lw lw*0.6 dt '-' lc rgb "#ff9980" title '',\
    [9e-6:r_c_neutrals] (analysis?1:1/0) * density_freefall(x, r_c_neutrals, rho_sonicpoint) axes x1y1 with lines lw lw*0.6 dt '-' lc rgb "#737373" title '',\
    [2.4e-4:1] (analysis?1:1/0) * density_hydrostatic(x, r_c_neutrals, rho_sonicpoint) axes x1y1 with lines lw lw*0.6 dt '-' lc rgb "#1a0cb5" title ''
    #[2.4e-4:1] (analysis?1:1/0) * density_updated_hydrostatic(x, r_c, rho_sonicpoint, rho_inf) axes x1y1 with lines lw lw*0.6 dt '-' lc rgb "#1a0cb5" title ''
   
}

system("mv *.png evolution")
#system("cd evolution/newplots; convert -delay 10 -loop 0 *.png evolution.gif")
