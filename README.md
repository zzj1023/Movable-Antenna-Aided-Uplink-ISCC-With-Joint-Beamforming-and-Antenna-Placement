# Movable-Antenna-Aided Uplink ISCC With Joint Beamforming and Antenna Placement

This repository contains the MATLAB implementation used to evaluate an
uplink integrated sensing, communication, and computation (ISCC) system
with movable antennas (MAs). The proposed method jointly optimizes MA
positions and transmit beamforming. Particle swarm optimization (PSO) is
used for the outer MA-position search, while SDR-SCA is used for the inner
beamforming problem.

The repository also includes the following benchmark schemes:

- Genetic algorithm (GA) based MA-position optimization
- Fixed-position antenna array (FPA)
- Greedy antenna selection (AS)
- Random-position antennas (RPA)

All schemes use the same channel realizations, sensing requirements,
computing parameters, and SDR-SCA beamforming solver.

## Requirements

- MATLAB
- CVX 2.2 or a compatible CVX release
- A CVX-supported SDP solver

After installing CVX, run `cvx_setup` and add this repository to the MATLAB
path:

```matlab
addpath(genpath(pwd));
```

## Repository Structure

```text
channel/            Fixed channel realizations for the main experiments
scenario_saved/     Nested fixed scenarios for different device counts
Plot/               Figure-generation scripts
res_P/              Results for the transmit-power sweep
res_M/              Results for the antenna-number sweep
res_Area/           Results for the MA-region-size sweep
res_K/              Results for the device-number sweep
sum_conve/          Convergence results
fig/                Exported figures
```

Core implementation files:

```text
run_SDR_sca_sum_rate_ISAC.m       SDR-SCA beamforming solver
channel_model_ISAC_deterministic.m
default_computing_config.m        Common computing and latency parameters
prepare_computing_constraints.m
optimize_ma_ga.m                  Common GA position optimizer
antenna_selection_greedy.m        Greedy AS implementation
```

## Reproducing the Results

Run MATLAB from the repository root. The scripts below regenerate the data
for each performance figure.

### Maximum transmit power

```matlab
ISCC_P
GA_P
FPA_P
AS_P
Random_P
```

### Number of receive antennas

```matlab
main_ISCC_M
GA_M
FPA_M
AS_M
Random_M
```

### MA receive-region size

```matlab
main_ISCC_Area
GA_Area
FPA_Area
AS_Area
Random_Area
```

### Number of devices

```matlab
main_ISCC_K
GA_K
FPA_K
AS_K
Random_K
```

### Convergence comparison

```matlab
main_ISCC
main_GA
main_FPA
main_AS
main_Random
```

The full PSO and GA experiments require many CVX solves and can take a long
time. For a quick code check, reduce the population size and number of
iterations temporarily, but restore the reported settings before producing
paper results.

## Plotting

After generating the result files, change to the `Plot` directory and run:

```matlab
plot_P
plot_M
plot_Area
plot_K
Convergence
```

## Reproducibility Notes

- Simulation scripts set explicit random seeds.
- The `channel` and `scenario_saved` directories contain fixed channel data.
- `default_computing_config.m` defines the common task sizes, CPU
  allocations, bandwidth, CPU budget, and latency threshold.
- Generated result MAT files are included so that figures can be reproduced
  without rerunning every optimization.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).

## Citation

If this code is useful in your research, please cite the associated paper.
The final bibliographic entry will be added after publication.
