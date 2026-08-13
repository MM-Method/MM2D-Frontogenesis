#!/bin/bash

# Exit script if an unhandled error occurs
set -e

# Store the absolute path to the executable and template
PROGRAM="/home/ahokpunn/multimoment/github/MM2D-main/MMM_1Dconvect"
TEMPLATE="$(pwd)/setup.templ"
COREperChiplet=8


# Ensure required files exist before starting
if [ ! -f "$TEMPLATE" ]; then
echo "Error: Template file '$TEMPLATE' not found!"
exit 1
fi

if [ ! -x "$PROGRAM" ]; then
echo "Error: Executable '$PROGRAM' not found or not executable!"
exit 1
fi

# Define schemes (Columns in Tab.A)
SCHEMES=("mmex4" "mmex6" "mmex8" "mmex10" "fdex4" "fdex6" "fdex8" "fdex10")

# Define grid points array indexed as [ACCURACY_INDEX, SCHEME_INDEX]
# Rows: 0=1e-2, 1=1e-3, 2=1e-4, 3=1e-5
GRID_POINTS=(
"160" "122" "114" "111" "380"  "268"  "244" "236"  # 1e-2
"302" "220" "186" "174" "728"  "478"  "398" "362"  # 1e-3
"550" "342" "272" "250" "1338" "750"  "584" "514"  # 1e-4
"970" "514" "378" "328" "2394" "1118" "816" "684"  # 1e-5
)

ACCURACIES=("1e-2" "1e-3" "1e-4" "1e-5")

# Benchmarking parameters: 16 parallel instances x 8 batches = 128 total runs per scheme
PARALLEL_PER_BATCH=16
NUM_BATCHES=8

# Pre-create output folders for each scheme
for scheme in "${SCHEMES[@]}"; do
mkdir -p "$scheme"
done

echo "=========================================================="
echo " Starting Benchmarks: Scheme-by-Scheme                    "
echo " (16 Parallel Instances x 8 Batches = 128 Runs Per Scheme) "
echo "=========================================================="

# 1. Outer Loop: Accuracy Levels (Runs sequentially)
for acc_idx in "${!ACCURACIES[@]}"; do
accuracy="${ACCURACIES[$acc_idx]}"

echo ""
echo "=========================================================="
echo " Starting Accuracy Level: $accuracy"
echo "=========================================================="

# 2. Middle Loop: Schemes (Runs sequentially, one scheme at a time)
for scheme_idx in "${!SCHEMES[@]}"; do
scheme="${SCHEMES[$scheme_idx]}"

# Calculate linear index for grid points array
flat_idx=$(( acc_idx * 8 + scheme_idx ))
grdpoint="${GRID_POINTS[$flat_idx]}"

# Determine CFL Rule
if [[ "$scheme" == fd* ]]; then
cfl="0.8"
else
cfl="0.4"
fi

# Determine NMETH Rule
if [[ "$scheme" == fd* ]]; then
nmeth="NC0"
else
if [[ "$accuracy" == "1e-5" ]]; then
if [[ "$scheme" == "mmex10" ]]; then
nmeth="NC6"
else
nmeth="NC4"
fi
else
nmeth="NC4"
fi
fi

# Generate local setup.txt inside the scheme's folder
sed \
-e "s/GRDPOINT/$grdpoint/g" \
-e "s/CFL/$cfl/g" \
-e "s/SCHEME/$scheme/g" \
-e "s/NMETH/$nmeth/g" \
"$TEMPLATE" > "$scheme/setup.txt"

echo "[LAUNCHING SCHEME] $scheme ($accuracy) -> Running 160 instances across $NUM_BATCHES batches..."

# 3. Inner Loop: Batches (1 to 16)
for batch in $(seq 1 $NUM_BATCHES); do

# Launch 8 parallel instances in the background for this batch
for inst in $(seq 1 $PARALLEL_PER_BATCH); do
# Unique run number from 1 to 128
run_num=$(( (batch - 1) * PARALLEL_PER_BATCH + inst ))

core=$((COREperChiplet * (inst-1)  ))

(
cd "$scheme"
log_file="run${run_num}_${accuracy}_grd${grdpoint}.txt"
taskset -c $core $PROGRAM  > $log_file 2>&1
) &  # <-- Spawns process in background
done

# Wait for all 16 background instances in this batch to complete
wait

done

echo "[FINISHED SCHEME] $scheme for $accuracy (128 output files generated)."

done

echo "----------------------------------------------------------"
echo " Completed all schemes for accuracy level: $accuracy"
echo "----------------------------------------------------------"

done

echo ""
echo "=========================================================="
echo " All Benchmarks Completed Successfully! "
echo "=========================================================="
