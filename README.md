# OrbitalMechanics

Software by Carlo Zambaldo (info@carlozambaldo.it) This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License.

I really want to thank also my project team: Andrea Zosi, Guglielmo Gomiero and Lorenzo Sekules

## ORBIT class
ORBIT is an object describing an orbit, the principal fields are:
    - a        : semi-major axis [km]
    - e        : eccentricity [-]
    - i        : inclination [rad]
    - O        : RAAN [rad]
    - w        : anomaly of pericentre [rad]
    - theta    : true anomaly (usually initial position) [rad]
    - mu       : gravitational constant, for Earth use astroConstants(13)
 
 in addition, the by defining the orbit using ORBIT(...) the following
 fields are filled in:
    - r0        : position on orbit at given theta
    - v0        : velocity on orbit at given theta
    - T         : period of the orbit [s]
    - energy    : energy of the orbit
    - p         : orbital parameter
    - h_vect    : angular momentum vector (per unit mass)
    - e_vect    : eccentricity vector
    - rA        : radius of apoapsis
    - rP        : radius of periapsis
    - theta_inf : true anomaly asymptote
    - delta     : deflection
    - v_inf     : hyperbolic eccess speed
    - type      : char [c,e,p,h] circular, elliptical, parabolic, hyperbolic
                 and [+,=,-] respectively prograde, polar, retrograde


### brief explanation of the class
 - DEFINITION OF AN ORBIT OBJECT
    the orbit can be defined in parametrical form by inserting all 7
    parameters like: obj = ORBIT(a,e,i,O,w,theta,mu)
    if the orbit is given in cartesian form use the function keplerian to
    define the orbit like: obj = ORBIT.keplerian(r0, v0, mu).
    To define the orbit in degrees use ORBIT.define_in_deg(a,e,i,O,w,theta,mu).
    WARNING: r_vect, v_vect and h are initialized as 0 vectors, if
    there is a problem in the script these values are likely to remain
    zero. Moreover, check the 'type' proprety, in these cases it usually
    displays an 'error' message.

 - DEFINITION OF A SIMPLIFIED ORBIT OBJECT
    it is also possible to only store a,e,i,O,w,theta,mu values by
    calling the obj = simpleOrbit(a,e,i,O,w,theta,mu) method. If this type
    of orbit has to be filled, in a second moment, with all the other
    parameters: use obj.fillOrbit.
    
 - PLOTTING OF AN ORBIT
    to plot a 3D orbit use the method obj.plotOrbit [if needed add the
    optional parameters], to propagate an orbit use instead
    obj.propagateOrbit, this last method uses ode113 to propagate the
    orbit given a tspan, it is possible to modify the propagation
    method (i.e. if perturbances have to be taken into account) using 
    option.propagationType. Write help propagateOrbit to know more.

 - PLOTTING OF GROUND TRACKS
    ground-tracks plotting can simply be performed by using
    computeGroundTrack method without requiring any output. If
    computeGroundTrack is called requiring the outputs use
    plotGroundTrack to also plot the required groundTrack.
     note: groundTrack function is deprecated and in a future release
     will no longer work.

 - CONVERT BETWEEN CARTESIAN AND KEPLERIAN ELEMENTS AND VICE VERSA
    two possible approaches are available: if the orbit is not
    hyperbolic nor parabolic the user is highly encouraged to use
    kep2car and car2kep methods. The second approach is
    object-oriented and requires the definition of an ORBIT object to
    work: use therefore cartesian or keplerian methods to convert
    between the two rapresentations.

 - EXTRAS
     to compute the time of flight between two points use TOF()
     to check if a position vector is on the orbit use isOnOrbit
 

### list of methods
 methods:
  - [t,y] = propagateOrbit(obj,tspan,options,grafica,r_segnato);
  - [r_prop_vect] = plotOrbit(obj,theta_vect,d_theta,grafica,theta_segnato);
  - [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)
  - type = typedef(obj)
  - delta_time = TOF(obj, thetas);
  - delta_theta = AOF(obj, delta_time, toll)
  - v_esc = escapeVelocity(obj, position)
  - check = isOnOrbit(obj, r_vect)
  - obj = fillOrbit(obj)

 static methods:
  - obj = define_in_deg(a,e,i_deg,O_deg,w_deg,theta_deg,mu)
  - obj = simpleOrbit(a,e,i,O,w,theta,mu)
  - [r_vect, v_vect]  = kep2car(a,e,i,O,w,theta,mu);
  - [a,e,i,O,w,theta] = car2kep(r_vect, v_vect, mu);
  - [r_vect, v_vect] = cartesian(obj, theta);
  - [obj] = keplerian(r_vect, v_vect, mu);
  - [] = plotGroundTrack(lon, lat, options)



> [!WARNING]
> THE ORBIT OBJECT TO WORK PROPERLY HAS TO BE IN A FOLDER NAMED
> PRECISELY "@ORBIT" WITH ALL ITS METHODS, THIS FOLDER HAS INDEED 
> TO BE ADDED TO THE PATH. NO EXCEPTIONS CAN BE MADE.
