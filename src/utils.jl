using CSV
using GeometryBasics
using LinearAlgebra
using GLMakie
using GLMakie.GLFW
using GLMakie: to_native
using FileIO
using ColorSchemes
using DataFrames
using OrdinaryDiffEq
using JLD2


GLMakie.activate!()

## Define constants:
#----------------------------------------------------------------------
	small             = 1.0e-10
	infinite          = 999999.9
	undefined         = 999999.1
	rad               = 180.0 / pi
	twopi             = 2.0 * pi
	halfpi            = pi * 0.5
	ft2m              = 0.3048
	mile2m            = 1609.344
	nm2m              = 1852
	mile2ft           = 5280
	mileph2kmph       = 0.44704
	nmph2kmph         = 0.5144444
	re                = 6378.137         # km
	flat              = 1.0/298.257223563
	omegaearth        = 7.292115e-5      # rad/s
	mu = MU = μ       = 398600.4418      # km3/s2
	mum               = 3.986004418e14   # m3/s2
	eccearth          = sqrt(2.0*flat - flat^2)
	eccearthsqrd      = eccearth^2
	renm              = re / nm2m
	reft              = re * 1000.0 / ft2m
	tusec             = sqrt(re^3/mu)
	tumin             = tusec / 60.0
	tuday             = tusec / 86400.0
	omegaearthradptu  = omegaearth * tusec
	omegaearthradpmin = omegaearth * 60.0
	velkmps           = sqrt(mu / re)
	velftps           = velkmps * 1000.0/ft2m
	velradpmin        = velkmps * 60.0/re
	degpsec           = (180.0 / pi) / tusec
	radpday           = 2.0 * pi * 1.002737909350795
	speedoflight      = 2.99792458e8    # m/s
	au                = 149597870.0    	# km
	earth2moon        = 384400.0 		# km
	moonradius        = 1738.0 		    # km
	sunradius         = 696000.0 		# km
	masssun           = 1.9891e30
	massearth         = 5.9742e24
	massmoon          = 7.3483e22
#----------------------------------------------------------------------


"""
    Orbit(a, e, i, Ω, ω, M0, ν, p, u, l, Π, n, P, μ, r0, v0, h0)

    Creates a orbit structure containing all COEs and extras
"""
struct Orbit
   a::Float64
   e::Float64
   i::Float64
   Ω::Float64
   ω::Float64
   M0::Float64
   ν::Float64
   p::Float64
   u::Float64
   l::Float64
   Π::Float64
   n::Float64
   P::Float64
   μ::Float64
   r0::Vector{Float64}
   v0::Vector{Float64}
   h0::Vector{Float64}
end

"""
    genSatellite(a, e, i, Ω, ω, M0)

    Function used to generate a single orbit structure given some COEs
"""
function genOrbit(a, e, i, Ω, ω; M0=nothing, ν=nothing)

   if M0 === nothing && ν === nothing
    error("Must define either M0 or ν (ν takes precedence if both defined)")
   end
   if ν === nothing
    ν = newtonm(e, M0)[2]
   else
    M0 = 0
   end

   p = a * (1 - e^2)
   u = ν + ω
   l = Ω + ω + ν
   Π = Ω + ω
   n = sqrt(μ / a^3)
   P = 2.0 * π / n

   # Define values here to make calculating r0 and v0 easier
   cfw = cos(u)
   sfw = sin(u)
   cOm = cos(Ω)
   sOm = sin(Ω)
   ci = cos(i)
   si = sin(i)
   cw = cos(ω)
   sw = sin(ω)

   # Initial position vector
   r0 = p/(1 + e*cos(ν)) * [cfw*cOm - ci*sfw*sOm; ci*cOm*sfw + cfw*sOm; si*sfw]

   # Initial velocity vector
   v0 = sqrt(μ/(a*(1-e^2))) *
      [-cOm*sfw - sOm*ci*cfw - e*(cOm*sw + sOm*cw*ci);
      cOm*ci*cfw - sOm*sfw - e*(sOm*sw - cOm*cw*ci);
      si*(cfw + e*cw)]
   # Initial angular momentum
   h0 = cross(r0, v0)

   # Create an Orbit structure with these elements and return it
   return Orbit(a, e, i, Ω, ω, M0, ν, p, u, l, Π, n, P, μ, r0, v0, h0)
end


function genOrbit(R, V)
    appropriate_sizes = [(3,), (1, 3), (3, 1)]
    if size(R) ∉ appropriate_sizes || size(V) ∉ appropriate_sizes
        error("R and V must contain three elements")
    end

    a, e, i, Ω, ω, ν = elOrb(R,V)
    return genOrbit(a, e, i, Ω, ω; ν=ν)

end

function elOrb(R,V)
    appropriate_sizes = [(3,), (1, 3), (3, 1)]
    if size(R) ∉ appropriate_sizes || size(V) ∉ appropriate_sizes
        error("R and V must contain three elements")
    end

    H, N, E = hnevec(R, V)
    i, Ω, ω, ν = angles(R, V, H, N, E)
    a, e = sizeshape(R, V, E)
    return a, e, i, Ω, ω, ν
end

function hnevec(R, V)
    appropriate_sizes = [(3,), (1, 3), (3, 1)]
    if size(R) ∉ appropriate_sizes || size(V) ∉ appropriate_sizes
        error("R and V must contain three elements")
    end
    
    H = cross(R, V)
    N = cross([0; 0; 1], h)
    E = (1/μ) * [(norm(V)^2 - (μ / norm(R))) * R - (dot(R, V)) * V]

    return H, N, E
end

function angles(R, V, H, N, E)
    i = vecangle([0; 0; 1], H)
    Ω = N[2] >= 0 ? vecangle([1; 0; 0], N) : 2π - vecangle([1; 0; 0], N)
    ω = E[3] >= 0 ? vecangle(N, E) : 2π - vecangle(N, E)
    ν = dot(R, V) >= 0 ? vecangle(E, R) : 2π - vecangle(E, R)
    return i, Ω, ω, ν
end

vecangle(v1, v2) = acos((dot(v1, v2)) / (sqrt(dot(v1, v1)) * sqrt(dot(v2, v2))))

function sizeshape(R, V, E)
    a = μ / (-norm(V)^2 + (2μ) / norm(R))
    ecc = norm(E)
    return a, ecc
end

function r̈(ẋ, x, g, t)
    x_norm = norm(x[1:3])
    ẋ[1] = x[4]
    ẋ[2] = x[5]
    ẋ[3] = x[6]
    ẋ[4] = -μ/x_norm^3 * x[1]
    ẋ[5] = -μ/x_norm^3 * x[2]
    ẋ[6] = -μ/x_norm^3 * x[3]
end

function hohmannTransfer(orbit1::Orbit, orbit2::Orbit)

    Rpt = orbit1.a * (1 - orbit1.e)
    Rat = norm(orbit2.r0)
    at = (Rpt + Rat) / 2
    # v0t = [sqrt(2*(μ/orbit_tran.r0[1] - μ/(2*orbit_tran.a)));
    #                  sqrt(2*(μ/orbit_tran.r0[2] - μ/(2*orbit_tran.a)));
    #                  sqrt(2*(μ/orbit_tran.r0[3] - μ/(2*orbit_tran.a)))]
    et = 1 - Rpt / at

    return genOrbit(at, et, orbit1.i, orbit1.Ω, orbit1.ω; ν=orbit1.ν)

end

function getOrbit(orbit::Orbit; numP=1, refOrb=nothing)
    tf = orbit.P + 1
    steps = 300
    if refOrb != nothing
        Δt = refOrb.P / steps
    else
        Δt = orbit.P / steps
    end

    X₀ = [orbit.r0..., orbit.v0...]
    
    ode = ODEProblem(r̈, X₀, (0.0, tf))
    sol = solve(ode, Tsit5(), reltol=1e-5, abstol=1e-5, saveat=Δt)
    u = hcat(sol.u...)'
    t = sol.t
    orb_whole = Int(floor(numP))
    orb_frac = numP % orb_whole
    r1 = []
    r2 = []
    r3 = []
    for i in 1:numP
        append!(r1, u[:,1])
        append!(r2, u[:,2])
        append!(r3, u[:,3])
        append!(t, i * t[end] .+ t)
    end
    return r1, r2, r3, t
end

function up_to_date(orbit1_params, orbit2_params, filename)
    if !isfile(filename)
        return false, nothing, nothing
    end
    orbit1 = load(filename, "orbit1")
    orbit2 = load(filename, "orbit2")
    numP1 = load(filename, "numP1")
    numP2 = load(filename, "numP2")
    if (orbit1.a != orbit1_params[1] || 
        orbit1.e != orbit1_params[2] ||
        orbit1.i != orbit1_params[3] ||
        orbit1.Ω != orbit1_params[4] ||
        orbit1.ω != orbit1_params[5] ||
        orbit1.ν != orbit1_params[6] ||
        numP1    != orbit1_params[7])

        return false, nothing, nothing

    elseif (orbit2.a != orbit2_params[1] || 
            orbit2.e != orbit2_params[2] ||
            orbit2.i != orbit2_params[3] ||
            orbit2.Ω != orbit2_params[4] ||
            orbit2.ω != orbit2_params[5] ||
            orbit2.ν != orbit2_params[6] ||
            numP2    != orbit2_params[7])

            return false, nothing, nothing
    end

    return true, orbit1, orbit2

end

"""
    newtonm(ecc, m)

    Calculates the eccentric and true anomalies using Newton-Raphson method.
    This function was created by Vallado
"""
function newtonm(ecc, m)
    numiter = 50
    small = 0.00000001
    halfpi = π * 0.5

    # Hyperbolic
    if ( (ecc-1.0 ) > small )
           # -------------------  initial guess -----------------------
            if ( ecc < 1.6  )
                if ( ((m<0.0 ) && (m>-pi)) || (m>pi) )
                    e0 = m - ecc
                  else
                    e0 = m + ecc
                end
              else
                if ( (ecc < 3.6 ) && (abs(m) > pi) )
                    e0 = m - sign(m)*ecc
                  else
                    e0 = m/(ecc-1.0 )
                end
            end
            ktr = 1
            e1 = e0 + ( (m-ecc*sinh(e0)+e0) / (ecc*cosh(e0) - 1.0 ) )
            while ((abs(e1-e0)>small ) && ( ktr<=numiter ))
                e0= e1
                e1= e0 + ( ( m - ecc*sinh(e0) + e0 ) / ( ecc*cosh(e0) - 1.0  ) )
                ktr = ktr + 1
            end
            # ----------------  find true anomaly  --------------------
            sinv= -( sqrt( ecc*ecc-1.0  ) * sinh(e1) ) / ( 1.0  - ecc*cosh(e1) )
            cosv= ( cosh(e1) - ecc ) / ( 1.0  - ecc*cosh(e1) )
            ν  = atan( sinv,cosv )
          else
            # --------------------- parabolic -------------------------
            if ( abs( ecc-1.0  ) < small )
                 #c = [ 1.0/3.0; 0.0; 1.0; -m]
                 #[r1r] = roots (c)
                 #e0 = r1r
                 s = 0.5  * (halfpi - atan( 1.5 * m ) )
                 w = atan( tan( s )^(1.0 /3.0 ) )
                 e0 = 2.0 * cot(2.0 * w)
                ktr = 1
                ν = 2.0  * atan(e0)
              else
                # -------------------- elliptical ----------------------
                if ( ecc > small )
                    # -----------  initial guess -------------
                    if ( ((m < 0.0 ) && (m > -pi)) || (m > pi) )
                        e0 = m - ecc
                      else
                        e0 = m + ecc
                    end
                    ktr = 1
                    e1 = e0 + (m - e0 + ecc * sin(e0)) / (1.0  - ecc*cos(e0))
                    while (( abs(e1-e0) > small ) && ( ktr <= numiter ))
                        ktr += 1
                        e0 = e1
                        e1 = e0 + (m - e0 + ecc*sin(e0)) / (1.0  - ecc*cos(e0))
                    end
                    # -------------  find true anomaly  ---------------
                    sinv = (sqrt(1.0 - ecc*ecc) * sin(e1)) / (1.0 - ecc*cos(e1))
                    cosv = (cos(e1) - ecc) / (1.0  - ecc*cos(e1))
                    ν  = atan(sinv, cosv)
                  else
                    # -------------------- circular -------------------
                    ktr = 0
                    ν = m
                    e0 = m
                end
            end
        end
    return e0, ν
end

"""
    wrapAngleRad(θ)

    Wraps input angle in radians to be between 0 and 2π.
"""
function wrapAngleRad(θ)

    if θ >= 2*π
        while θ >= 2*π
            θ -= 2*π
        end
    elseif θ < 0
        while θ < 0
            θ += 2*π
        end
	end

    return θ
end

"""
    wrapAngleDeg(θ)

    Wraps input angle in degrees to be between 0 and 360.
"""
function wrapAngleDeg(θ)

    if θ >= 360
        while θ >= 360
            θ -= 360
        end
    elseif θ < 0
        while θ < 0
            θ += 360
        end
	end

    return θ
end

