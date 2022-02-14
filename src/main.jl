#####################################################################################################
#
# Hohmann Transfer and Ground Track Visualization
#
# Description: This program takes two user-defined orbits and calculates the correct
# Hohmann transfer orbit.  It plots the two orbits and the transfer orbit
# along with a 3D model of Earth.  The program also shows how the 
# spacecraft progresses through the orbit as well as how the ground track
# extends across the surface of the Earth.
#
#####################################################################################################

now = time()

#####################################################################################################
# Define Orbits
#####################################################################################################
# a    = semi-major axis (km)
# e    = eccentricity (unitless)
# i    = inclination (radians)
# Ω    = right ascension of the ascending node (radians)
# ω    = argument of perigee (radians)
# ν    = true anomaly (radians)
# numP = number of periods to show in simulation

# Define the desired starting orbit
a1 = 16000; e1 = 0.0; i1 = deg2rad(10); Ω1 = deg2rad(0); ω1 = deg2rad(180); ν1 = deg2rad(0); numP1 = 1
# Define the desired target orbit
a2 = 9000; e2 = 0.1; i2 = deg2rad(10); Ω2 = deg2rad(0); ω2 = deg2rad(0); ν2 = deg2rad(0); numP2 = 2

save_video = true

#####################################################################################################
# Setup
#####################################################################################################
force_reload_all = false
force_reload_mesh = false
force_reload_orbits = false

println("Initializing                                             [0.000 sec]")
# Import the utility functions from utils.jl
include("utils.jl")

# Set the theme for Makie, these are most of the possible arguments
set_theme!(
	font = "Dejavu Serif",
	fonstize = 14,
	backgroundcolor = :black,
    color = :black,
	colormap = :viridis,
	markersize = 0.1,
    linestyle = nothing,
    visible = true,
    clear = true,
    show_axis = false,
    show_legend = false,
    scale_plot = true,
    center = true,
	update_limits = true,
	rowgap = 10,
    colgap = 20,
)

# Get a color palette, so we can use the colors later
colors = ColorSchemes.:seaborn_muted

# Load the Earth 3D model and texture if not already loaded
println("    Loading Earth 3D model and texture                   [$(round(time() - now, digits=3)) sec]")
if force_reload_all || force_reload_mesh || !@isdefined earth
    earth = load("../assets/globe.obj")
    earth_txt = load("../assets/8k_earth_daymap.jpg")
end

# Load the skybox 3D model and texture if not already loaded
println("    Loading Skybox 3D model and texture                  [$(round(time() - now, digits=3)) sec]")
if force_reload_all || force_reload_mesh || !@isdefined skybox
    skybox = load("../assets/sphere.obj")
    skybox_txt = load("../assets/8k_stars_milky_way.jpg")
end

# Load the cassini 3D model if not already loaded
println("    Loading spacecraft 3D model                          [$(round(time() - now, digits=3)) sec]")

if force_reload_all || force_reload_mesh || !@isdefined cassini
    cassini = load("../assets/cassini_small.stl")
end

# Initialize the figure with 4k size
fig = Figure()

# Create an observable that will keep track of the current time of the scene
scene_time = Observable(0.0)
# Create a function that calculates the current time using the scene_time observable
getTime = @lift(string("t = ", Int(round(($scene_time))), " sec;  ", round(($scene_time) / 3600, digits=2), " hr"))
# Create the location on the window that will display the current time
time_display = fig[1, 1] = Label(fig, getTime, textsize=40, color=:white)


# Create a button grid for all of the simulation speed buttons
fig[2, 1] = grid2 = GridLayout(tellwidth = false)
grid2[1, 1] = buttongrid = GridLayout(tellwidth = false)

# Define what values the buttons should have
button_values = (1, 2, 5, 7, 10, 15, 20)
# Assign the values to these buttons
buttons = b1, b2, b3, b4, b5, b6, b7 =
    buttongrid[1, 1:length(button_values)] =
    [Button(fig, label = "×$i") for i in button_values]

# Create a step size observable
step_size = Observable(1)

# Update the step size observable when a button is clicked
for (ix, b) in enumerate(buttons) 
    on(b.clicks) do val 
        step_size[] = button_values[ix]       
    end 
end


# Display a text box showing the currently selected speed
speed = @lift(string("Speed: ×", Int(round($step_size))))
Speed = buttongrid[1, 8] = Label(fig, speed, textsize=20, color=:white)

grid2[1, 2] = controlgrid = GridLayout(tellwidth = false)
# Create an observable specifying whether to rewind the visualization
rewind = Observable(false)
# Define what to display on the button based on the rewind observable
rewindScene = @lift($rewind ? "⏩" : "⏪")
# Create the rewind button
rewindButton = buttongrid[1, 9] = Button(fig, label=rewindScene)
# Update the button when it is clicked
on(rewindButton.clicks) do clicks
    rewind[] = !rewind[]
end


# Create an observable that specifies whether to pause the visualization
pause = Observable(false)

# Define the value to display on the pause button based on the pause observable
paused = @lift($pause ? "▷" : "▮▮")

# Create the pause button
pauseButton = buttongrid[1, 10] = Button(fig, label=paused)
# Update the pause button when it is clicked
on(pauseButton.clicks) do clicks
    pause[] = !pause[]
end

# Create an observable that specifies whether to pause the visualization
show_track = Observable(true)

# Define the value to display on the show track button based on the show_track observable
showing_track = @lift($show_track ? "showing track" : "hiding track")

# Create the show track button
showButton = buttongrid[1, 11] = Button(fig, label=showing_track)
# Update the show track button when it is clicked
on(showButton.clicks) do clicks
    show_track[] = !show_track[]
end

colsize!(grid2, 2, Relative(1/10))
# Initialize another scene, which will contain the 3D plot
scene3D = fig[3:10, 1] = LScene(fig, show_axis=false)
cam3d_cad!(scene3D.scene, near=1f-5)
colsize!(fig.layout, 1, Relative(1))


#####################################################################################################
# Rotation
#####################################################################################################
# Angular velocity of the Earth in rad/s2
ωₑ = 7.292115e-5
# Euler axis about which the Earth rotates (z-axis)
e = [0; 0; 1]

# Function that takes eigenaxis and angle and converts to quaternion
eigenaxisToQuat(e, α) = Quaternionf(e[1]*sin(α/2), e[2]*sin(α/2), e[3]*sin(α/2), cos(α/2))
# Function that calculates a skew symmetric matrix given a vector
skew(x) = [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]
# Function that takes an eigenaxis and angle and converts to a DCM
eigenaxisToDCM(e, α) = cos(α)I + (1 - cos(α)) * e * transpose(e) - sin(α) * skew(e)

# Function that takes a DCM and converts it to a quaternion
function DCMtoQuat(C)
    q = zeros(4)
    q[4] = 0.5 * (1 + C[1, 1] + C[2, 2] + C[3, 3]) ^ 0.5
    q[1:3] = 1 / (4 * q[4]) .* [C[2, 3] - C[3, 2]; C[3, 1] - C[1, 3]; C[1, 2] - C[2, 1]]
    return Quaternionf(q[1], q[2], q[3], q[4])
end

# Functions that calculate the elementary rotation DCMs (about x, y, z)
C₁(θ) = [1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
C₂(θ) = [cos(θ) 0 -sin(θ); 0 1 0; sin(θ) 0 cos(θ)]
C₃(θ) = [cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]

# Constant matrices that represent 90 deg rotations about x, y, and z axes
DCMx90 = C₁(π/2)
DCMy90 = C₂(π/2)
DCMz90 = C₃(π/2)


#####################################################################################################
# Generate, Save, or Load orbit and ground track data
#####################################################################################################
# Name of the file that will contain orbital data
filename = "datafile.jld"

# Store the orbital data into lists
orbit1_params = [a1, e1, i1, Ω1, ω1, ν1, numP1]
orbit2_params = [a2, e2, i2, Ω2, ω2, ν2, numP2]

# Define the radius of the Earth
Re = 6378.137

# Check if the requested orbits are valid
if a1 < Re || a2 < Re
    throw("Semi-Major axis cannot be smaller than the radius of the Earth")
elseif 1 < e1 < 0 || 1 < e2 < 0
    throw("Semi-Major axis must be between 0 and 1")
end

# Check whether new data needs to be generated, or if old data can be used
println("Checking if data is up-to-date                           [$(round(time() - now, digits=3)) sec]")

updated, orbit1, orbit2 = up_to_date(orbit1_params, orbit2_params, filename)

# If new data needs to be generated:
if !updated || force_reload_all || force_reload_orbits
    println("    Detected different orbits, generating new data       [$(round(time() - now, digits=3)) sec]")

    # Generate the orbital data based on the desired orbital parameters
    orbit1 = genOrbit(a1, e1, i1, Ω1, ω1; ν=ν1)
    orbit2 = genOrbit(a2, e2, i2, Ω2, ω2; ν=ν2)
    # Generate the Hohmann transfer orbit based on the two desired orbits
    orbit_tran = hohmannTransfer(orbit1, orbit2)

    # Extract the orbital data from the generated orbits
    r1₁, r2₁, r3₁, t₁ = getOrbit(orbit1; numP=numP1)
    r1₂, r2₂, r3₂, t₂ = getOrbit(orbit2; numP=numP2, refOrb=orbit1)
    r1ₜ, r2ₜ, r3ₜ, tₜ = getOrbit(orbit_tran, refOrb=orbit1)

    # Scale down the data by the radius of the Earth to match the size of the Earth 3D model
    scale = Re

    r1₁  = r1₁ ./ scale
    r2₁  = r2₁ ./ scale
    r3₁  = r3₁ ./ scale

    r1₂  = r1₂ ./ scale
    r2₂  = r2₂ ./ scale
    r3₂  = r3₂ ./ scale

    # Calculate the number of steps for half an orbit
    n = Int(round(maximum(size(r1ₜ)) / 2))

    # For the transfer orbit, only use the first half of the data and scale it down
    r1ₜ  = r1ₜ[1:n] ./ scale
    r2ₜ  = r2ₜ[1:n] ./ scale
    r3ₜ  = r3ₜ[1:n] ./ scale
    tₜ = tₜ[1:n]

    # Combine the starting orbit, the transfer orbit, and the target orbit into one vector for each axis
    r1_cmb = [r1₁..., r1ₜ..., r1₂...]
    r2_cmb = [r2₁..., r2ₜ..., r2₂...]
    r3_cmb = [r3₁..., r3ₜ..., r3₂...]

    # Combine the time vector to match the combined position vectors
    t_cmb = [t₁..., (tₜ .+ t₁[end])..., (t₂ .+ (tₜ[end] + t₁[end]))...]

    # Calculate the magnitude of the combined position vector at each point in time
    r_cmb_mag = [norm([r1_cmb[idx], r2_cmb[idx], r3_cmb[idx]]) for idx in 1:maximum(size(r1_cmb))]

    # Check if the requested orbit is valid
    if any(i -> i < 1, r_cmb_mag)
        println("Spacecraft reaches ", round(minimum(r_cmb_mag) * Re, digits=2), " km at its closest point")
        throw("Invalid orbit requested. Orbit intersects the Earth")
    end

    # Generate the unit position vectors for the combined orbit
    r1_cmb_unit = r1_cmb ./ r_cmb_mag
    r2_cmb_unit = r2_cmb ./ r_cmb_mag
    r3_cmb_unit = r3_cmb ./ r_cmb_mag

    # Initialize the ground track vectors to zero
    gnd_trk_x = zeros(size(r1_cmb_unit))
    gnd_trk_y = zeros(size(r2_cmb_unit))
    gnd_trk_z = zeros(size(r3_cmb_unit))

    # Calculate the difference between time steps
    Δt = t_cmb[2] - t_cmb[1]

    # Transform the unit position vector at each time step so that it "sticks to the Earth"
    for (ix, (t, r1, r2, r3)) in enumerate(zip(t_cmb, r1_cmb_unit, r2_cmb_unit, r3_cmb_unit))
        # Update the most recent position in the ground track with the unit position vector of the spacecraft
        gnd_trk_x[ix] = r1
        gnd_trk_y[ix] = r2
        gnd_trk_z[ix] = r3
        # Rotate the ground track by the rotation of the Earth in the given time span
        local gnd_trk = [gnd_trk_x[1:ix] gnd_trk_y[1:ix] gnd_trk_z[1:ix]] * eigenaxisToDCM(e, Δt * ωₑ)
        # Decompose and save the newly transformed ground track into each axis
        gnd_trk_x[1:ix] = [gnd_trk[1:ix, 1]...]
        gnd_trk_y[1:ix] = [gnd_trk[1:ix, 2]...]
        gnd_trk_z[1:ix] = [gnd_trk[1:ix, 3]...]
        
    end

    # Calculate the offset angle, θ, between the spacecraft and the ground track
    vec1 = [gnd_trk_x[1] gnd_trk_y[1]]
    vec2 = [r1_cmb_unit[1] r2_cmb_unit[1]]
    θ = atan(vec2[2], vec2[1]) - atan(vec1[2], vec1[1])

    # Undo the offset by rotating the ground track by the angle θ
    gnd_trk = [gnd_trk_x gnd_trk_y gnd_trk_z] * C₃(θ)
    # Save the newly transformed ground track
    gnd_trk_x = [gnd_trk[1:end, 1]...]
    gnd_trk_y = [gnd_trk[1:end, 2]...]
    gnd_trk_z = [gnd_trk[1:end, 3]...]
         
    # Save each relevant variable into a single JLD datafile
    jldsave(filename; orbit1, orbit2, numP1, numP2,
         r1₁, r2₁, r3₁, t₁,
         r1₂, r2₂, r3₂, t₂,
         r1ₜ, r2ₜ, r3ₜ, tₜ,
         r1_cmb, r2_cmb, r3_cmb, t_cmb, r_cmb_mag,
         r1_cmb_unit, r2_cmb_unit, r3_cmb_unit,
         gnd_trk_x, gnd_trk_y, gnd_trk_z)

    println("    New data successfully generated and saved            [$(round(time() - now, digits=3)) sec]")

else
    # If there is no need to generate new data, use the old data stored in the datafile
    println("    No change to orbits, using previously generated data [$(round(time() - now, digits=3)) sec]")

    # Load each variable from the datafile if not already loaded
    if !@isdefined r1₁

        orbit1, orbit2, numP1, numP2,
        r1₁, r2₁, r3₁, t₁,
        r1₂, r2₂, r3₂, t₂, 
        r1ₜ, r2ₜ, r3ₜ, tₜ, 
        r1_cmb, r2_cmb, r3_cmb, t_cmb, r_cmb_mag,
        r1_cmb_unit, r2_cmb_unit, r3_cmb_unit, 
        gnd_trk_x, gnd_trk_y, gnd_trk_z = 
            load(filename,"orbit1","orbit2","numP1","numP2",
            "r1₁","r2₁","r3₁","t₁",
            "r1₂","r2₂","r3₂","t₂",
            "r1ₜ","r2ₜ","r3ₜ","tₜ",
            "r1_cmb","r2_cmb","r3_cmb","t_cmb","r_cmb_mag",
            "r1_cmb_unit","r2_cmb_unit","r3_cmb_unit",
            "gnd_trk_x","gnd_trk_y","gnd_trk_z")
    end

    println("    Old data successfully loaded                         [$(round(time() - now, digits=3)) sec]")
end

#####################################################################################################
# Plot the data
#####################################################################################################

# Define the change in time and change in angle for each timestep
Δt = t₁[2] - t₁[1]
Δθ = ωₑ * Δt

# Define a function for retrieving the ground track at each timestep
getGroundTrack(t, length) = return [gnd_trk_x[1:t] gnd_trk_y[1:t] gnd_trk_z[1:t]]

# Create another observable for keeping track of the current step (integer)
step = Observable(0)

# Display the scene
glfw_window = to_native(display(fig))

# Plot the Skybox mesh
skybox_mesh = mesh!(
    scene3D,
    skybox,
    color = skybox_txt
)

# Plot the Earth mesh
earth_plot = mesh!(
    scene3D,
    earth,
    color = earth_txt
)

# Rotate the earth mesh so that it is upright
GLMakie.rotate!(earth_plot, Quaternionf(0.7071, 0, 0, 0.7071))

# Plot the entire starting orbit
first_orbit = lines!(
    scene3D,
    r1₁, r2₁, r3₁,
    linewidth = 3,
    color = colors[3]
)

# Plot the entire target orbit
second_orbit = lines!(
    scene3D,
    r1₂, r2₂, r3₂,
    linewidth = 3,
    color = colors[2]
)

# Plot the entire transfer orbit
transfer_orbit = lines!(
    scene3D,
    r1ₜ, r2ₜ, r3ₜ,
    linewidth = 3,
    linestyle = :dot,
    color = colors[5]
)

# Plot the ground track so far
ground_track = lines!(
    scene3D,
    @lift(getGroundTrack($step, 1000)),
    linewidth = 3,
    color = colors[1],
    visible = @lift(Bool($show_track * ($step >= 1 ? ($step / $step) : 0)))
)

# Plot the spacecraft mesh
spacecraft = mesh!(
    scene3D,
    cassini,
    color = colors[4]
)

# Make the spacecraft 4 percent of original size so it doesn't look ridiculous
spacecraft_scale = 0.04 * ones(3)
GLMakie.scale!(spacecraft, spacecraft_scale...)

max_r = maximum(abs.(r_cmb_mag))
# Update the camera to view the scene from a reasonable distance
update_cam!(scene3D.scene, cameracontrols(scene3D.scene), Vec3f0(1.2 *max_r, 1.2 *max_r, 1.2 *max_r), Vec3f0(0, 0, 0))


#####################################################################################################
# Run the Simulation
#####################################################################################################
println("Running Simulation                                       [$(round(time() - now, digits=3)) sec]")

function simulate(i)

    while pause[]
        # If the window is closed, stop the simulation
        if !events(scene3D.scene).window_open[]
            break
        end

        # If the "Esc" button is pressed, close the window and stop the simulation
        if ispressed(fig.scene, Keyboard.escape) # Press Esc to close the window
            GLFW.SetWindowShouldClose(glfw_window, true) # this will close the window after all callbacks are finished
        end
        sleep(1/20)
    end

    if !rewind[]
        # Increment the current step by the step size
        step[] += step_size[]
        # If the simulation reaches max time, reset the step number to 1
        if step[] >= maximum(size(r1_cmb)) - step_size[]
            step[] = 1
            #continue
        end
    else
        # Decrement the current step by the step size
        step[] -= step_size[]
        # If the simulation reaches zero time, reset the step number to max time
        if step[] <= 0 + step_size[]
            step[] = maximum(size(r1_cmb))
            #continue
        end
    end
    
    
    # Move the spacecraft mesh to the current position on the orbit
    GLMakie.translate!(spacecraft, r1_cmb[step[]], r2_cmb[step[]], r3_cmb[step[]])

    # Update the current scene time
    scene_time[] = step[] * Δt

    # Rotate the earth by the amount it should rotate in the given time
    GLMakie.rotate!(earth_plot, DCMtoQuat(DCMx90 * eigenaxisToDCM(e,  scene_time[] * ωₑ)))
    # Rotate the ground track to match the rotation of the Earth
    GLMakie.rotate!(ground_track, DCMtoQuat(eigenaxisToDCM(e,  scene_time[] * ωₑ)))
    # if !show_track[]
    #     local ground_track = 0 
    # end
    # Wait for 1/45 of a second to give the simulation time to load
    sleep(1/45)

    # If the "Esc" button is pressed, close the window and stop the simulation
    if ispressed(fig.scene, Keyboard.escape) # Press Esc to close the window
        GLFW.SetWindowShouldClose(glfw_window, true) # this will close the window after all callbacks are finished
    end

end

if save_video
    framerate = 30
    timestamps = range(0, maximum(size(r1_cmb)), step=step_size[])
    record(fig, "../hohmann-transfer.mp4", timestamps; framerate=framerate) do i
        simulate(i)
    end
else
    for orbit ∈ 0:10
        # Go through the entire combined orbit
        for i ∈ 1:step_size[]:maximum(size(r1_cmb))
            simulate(i)
            # If the window is closed, stop the simulation
            if !events(scene3D.scene).window_open[]
                break
            end
        end
        if !events(scene3D.scene).window_open[]
            break
        end
    end
end

# The simulation finishes, close the window
if events(scene3D.scene).window_open[]
    GLFW.SetWindowShouldClose(glfw_window, true)
end

println("    Simulation has ended                                 [$(round(time() - now, digits=3)) sec]")
