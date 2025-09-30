# PTRAC-to-PoliMi

A Python algorithm that uses the `h5py` package to read MCNP6.3's HDF5-formatted PTRAC output and convert it into the space-delimited collision file format required by the MCNPX-PoliMi post-processor (MPPost).

***

### Motivation

The MCNPX-PoliMi code, a modified version of MCNPX, stores detailed particle collision data in a specific file format for post-processing. This output is crucial for MPPost, a program that simulates detector responses and performs advanced data analysis. However, newer versions of MCNP, such as MCNP6.3, use a different format (HDF5) for their Particle Track (PTRAC) output. The `PTRAC-to-PoliMi` script was developed to create a compatible collision data file from this modern PTRAC output, allowing users to leverage MPPost with MCNP6.3 simulation results.

***

### Methodology

The script uses the `h5py` Python package to parse the MCNP6.3 HDF5 PTRAC output, focusing on the `ptrack` group, which contains datasets for all events in a history.

#### Data Extraction and Conversion
The conversion process involves these steps:

1.  **Event Ordering**: The script uses the `RecordLog` dataset to process events in the same order that MCNP simulates them. This chronological ordering is essential for accurately calculating energy deposition.
2.  **Energy Deposition**: The PTRAC `Collision` dataset provides the particle's energy *after* a collision, not the energy deposited. The script calculates the energy deposition for each collision event by taking the difference between the particle's energy before the event and its energy after the event.
3.  **Field Mapping**: The script maps relevant fields from the PTRAC datasets to the required PoliMi format. It uses data from the `Collision` dataset to populate fields such as `History_number`, `Reaction_Type`, `Target_nucleus`, `Time`, and `Weight`.


#### Handling Missing Data
While the MCNP6.3 PTRAC output is quite detailed, some fields required by the MCNPX-PoliMi format are handled differently or are not directly available.
* **"Generation number"**: This field is not available in the MCNP6.3 PTRAC output.
* **"Particle number"**: A custom particle tracking algorithm is required to fully replicate this. The script uses a simplified approach.
* **"Number scatterings"**: The `num_collisions_this_branch` field from the PTRAC `Collision` dataset is used as a replacement for this field.

***

### Verification and Results

The functionality of the PTRAC-to-PoliMi script has been verified by comparing its output with the native output from MCNPX-PoliMi across various test cases.

* **Energy Deposited**: Plots comparing the energy deposited show that the results from the PTRAC-derived data are highly consistent with the original PoliMi output. This confirms the accuracy of the energy deposition calculation method used in the script.

* **Particle Number and Multiplicity**: Comparisons of detection rates versus particle number also demonstrate strong agreement between the two methods, indicating the script's capability to preserve the overall particle multiplicity distribution.

* **Correlations**: Plots of time difference (cross-correlation) between detectors confirm that the PTRAC-to-PoliMi conversion accurately preserves the temporal correlations of particle events. This is a critical feature for advanced time-correlated measurements.


These results validate that the `PTRAC-to-PoliMi` script successfully generates a compatible and accurate input file for post-processing with MPPost.

***

### Usage

This script requires the `h5py` Python library. You can install it using `pip install h5py`.

To run the conversion, use the following Python function call:

```python
process_ptrac_hdf5_to_polimi_ascii('path/to/your/ptrac_output.h5', 'path/to/your/polimi_collision_data.txt')
