# R3B_sim_proton_tracking
Simple code to analyze output of the R3B simulation with FOOT and CALIFA detectors

# Download
git clone https://github.com/aldros/R3B_sim_proton_tracking.git

# Install
#### First source your R3BRoot install
cd R3BRoot_build
source config.sh
#### Then
cd R3B_sim_proton_tracking
make

# How to use
cd R3B_sim_proton_tracking
./run_proton_tracking --file=sim.root
"sim.root" is the output  of the R3BRoot simulation that should be located in input/
Result will be located in output/

