#!/usr/bin/env bash

# Default binary locations; override with CISTEM_BIN as needed.
path_to_cistem="${CISTEM_BIN:-/cluster/apps/cisTEM/2.0.0-alpha}"

# Use Singularity to invoke cisTEM inside the containerized environment.
function_to_execute="simulate"

# Note: Pass --nv when the host has GPUs available; combine with --bind for additional filesystems.

################################
# File config
################################

# Example input: ../2w0o_complex_ATOM.pdb

input_pdb_file=replace_here
# The parser tolerates standard PDB files but may fail on uncommon records.

input_basename=$(basename "${input_pdb_file}")
input_stem="${input_basename%.*}"
output_filename="${input_stem}__noise_${1}x${2}e_${3}t_${4}p_${6}d_${5}deg.mrc"

output_size=-500
# Negative values represent box size; larger boxes reduce cisTEM edge artifacts.

n_threads=16 # 3D simulations remain fast; thread count mainly affects IO waits.


################################
# Blur configuration
################################
linear_scaling_of_PDB_bfactors=1.0 # Ensure source BFactors remain realistic.
per_atom_bfactor=4.0 # Global BFactor boost to mimic motion blur.

################################
# Imaging parameters
################################
pixel_size=1.0
CS=2.7 # mm
KV=300 # kev
OA=100 # micron, objective aperture won't affect a 3D sim

wanted_defocus=$6 # Angstrom
minimum_thickness=$3

exposure_per_frame=$1 # e-/Ang^2; combined with frame count for total dose.
exposure_rate=3.0 # e-/pixel/s
n_frames=$2
pre_exposure=0 # Adjust if early frames were discarded experimentally.
n_particles=1
extraphase=$4 # 0.5 = pi/2
water_scaling=0.0 # Increase to add structured water noise.
angle=$5


################################
# Fixed arguments
################################
from_future_args=" --only-modify-signal-3d=2 --wgt=0.225 --water-shell-only"

################################
# Optional arguments
################################
# --save-detector-wavefunction : persist detector wavefunction probabilities.
# --skip-random-angles : reuse orientation for consistent comparisons.
# --max_number_of_noise_particles=6 : increase noise particle count.

optional_args="   "

################################
# Run the thing
################################

${path_to_cistem}/${function_to_execute}  --emulate_tilt_angle ${angle}  ${from_future_args} ${optional_args} << EOF
${output_filename}
$output_size
$n_threads
${input_pdb_file}
no
$pixel_size
$linear_scaling_of_PDB_bfactors
$per_atom_bfactor
no
$n_particles
$wanted_defocus
${extraphase}
${exposure_per_frame}
$exposure_rate
${n_frames}
yes
yes
../stars/tilt_${angle}.star
1
$water_scaling
0
$KV
$OA
$CS
0.0
$pre_exposure
32
$minimum_thickness
5
0.0
2.0
0.1
0.0001
0.0
0.0
0.0
0.0
EOF


