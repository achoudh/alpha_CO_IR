
# ### introduction ###
# # CO layer parameters
degrees::Float64  = pi/180

# Number of molecules in an unit cell of monolayer and overlayer
nmols_uc::Int64   = 4
nmols_ucol::Int64 = 8

# nmols_ml::Int64 = nmols_uc*nx*ny
# nmols_ol::Int64 = nmols_ucol*nx*ny*nz
# ndofs_ml::Int64 = 1 + 5*nmols_ml
# nmols_ol2::Int64 = nmols_ucol*nx*ny # 2 lowest layers of overlayer


r_CO::Float64 = 1.14e-10       # CO bondlength in m
a0_CO::Float64 = 5.64e-10      # CO layer lattice constant, m
a0_NaCl::Float64 = 5.64e-10    # NaCl lattice constant, m
a0_surf::Float64 = 3.99e-10    # NaCl surface lattice constant

v::Float64 = 0.4903e-10 #Oxygen
w::Float64 = 0.6437e-10 #Carbon
# include("C:/Users/achoudh/ownCloud/my work/CO_NaCl-estat/Estat_results/CO_NaCl_Estat/constants.jl")
include("./lattice_construction.jl")

z_ml = 3.35e-10/a0_surf

# construct an overlayer

# orientation of molecules in an overlayer's unit cell
θ_uc::Vector{Float64} = [ 3, 1, 3, 1, 1, 3, 1, 3]*pi/4.0  # The old structure [ 3, 1, 3, 1, 3, 1, 3, 1]*pi/4.0 is not correct
ϕ_uc::Vector{Float64} = [-1, 0,-1, 0, 2, 1, 2, 1]*pi/2.0 
trig_uc = (sin.(θ_uc), cos.(θ_uc), sin.(ϕ_uc), cos.(ϕ_uc))
# overlayer-surface distance (reduced units)
z_ol = z_ml + 0.5*a0_CO/a0_surf #+ 10.00*a0_CO/a0_surf
# get an overlayer molecules' reduced positions and orientation
com0_ol, phi_ol, theta_ol = overlayer(θ_uc, ϕ_uc, z_ol)
# δr_ol = zeros(Float64,nmols_ol2,2)


######################
# Exciton parameters #
######################

# Data to built up the wavenumber array
range   = 10  # "Wavenumbers"
dtponts = 5*200
step    = 2 * (range / dtponts)
# νk = collect(ν0[2] - range :step:ν0[1] + range)
Δν = 0.2 # cm-1 FWHM of the Gaussian convolution

μ00, μ11, μ01 = -0.112, -0.087, 0.10 # "Debyes"; μ00 and μ11: R.Disselkamp et al., Surface Science 240 (1990) 193-210; for 12C16O. μ01 calculated for 13C18O.
unit1         = 5034.12*1e-30 # conversion factor from Debye^2/m^3 to wavenumber
unit2         = 7.51691023

#electric field
θe       = 45.0*pi/180 # Tilt of incident beam
nar, ncr = 1.0, 1.52
nrat     = nar/ncr

Tx = 2*cos(θe)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2) / ((ncr/nar) + (cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));
Ty = 2 / (1 + (ncr/nar)*(cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));
Tz = 2*(sin(θe))^2 / (1 + (nar/ncr)*(cos(θe))^(-1)*sqrt(1 - ((nar/ncr)^2)*(sin(θe))^2));

Tp, Ts = [[sqrt(Tx), 0.0, sqrt(Tz)],[0.0, sqrt(Ty), 0.0]]

ep, es = Tp, Ts;