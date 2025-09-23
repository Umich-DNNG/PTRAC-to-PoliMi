import h5py
import numpy as np
from tqdm import tqdm
from sys import stdout
import concurrent.futures
import os

# This dictionary maps MCNP reaction type (MT) numbers to PoliMi NTYN codes.
# It's used to standardize reaction identifiers in the output.
mt2ntyn = { -1:1  ,
            -2:2  ,
            -4:5  ,
             2:-99,
           101:0  ,
           102:0  ,
           107:0   }
for i in range(51,92): mt2ntyn[i] = -1

# This dictionary stores the mass of specific nuclides (in AMU).
# It's used for energy deposition calculations in inelastic scattering.
masses = {1001:0.99916733,
          6012:11.9078563}

# This dictionary stores Q-values (in MeV) for specific reactions.
# It's used for energy deposition calculations, particularly for reactions like (n, alpha).
qvals = {6012:
           {107:-5.70205}}

def calculate_energy_deposition(prev_erg,prev_vec,ipt,rxn,za,erg,vec,x,y,z,tme,wgt,nps):
    """
    Calculates the energy deposited in a collision event based on particle type and reaction.
    
    Args:
        prev_erg (float): Energy of the particle before the collision.
        prev_vec (np.array): Direction vector of the particle before the collision.
        ipt (int): Particle type (1 for neutron, 2 for photon).
        rxn (int): Reaction type code.
        za (int): ZAID of the target nuclide.
        erg (float): Energy of the particle after the collision.
        vec (np.array): Direction vector of the particle after the collision.
        x, y, z, tme, wgt, nps: Other event data for context in warning messages.

    Returns:
        float: The calculated energy deposited in the collision.
    """
    # --- Neutron collisions ---
    if ipt == 1:
        # Elastic scatter: energy difference
        if rxn == 2:
            dep = prev_erg - erg
        # Inelastic scatter: calculated from kinematics
        elif (51 <= rxn) and (rxn <=91):
            d1 = prev_vec
            d2 = vec
            mu = np.dot(d1,d2) / (np.dot(d1,d1)*np.dot(d2,d2))**0.5
            try:
                dep = (prev_erg + erg - 2*mu*(prev_erg*erg)**0.5)/masses[za]
            except KeyError:
                stdout.write(f"WARNING: Inelastic Collision Missing Mass, {nps} {za}\n")
                dep = prev_erg
        # Radiative capture
        elif (rxn == 102) | (rxn == 101):
            try:
                dep = prev_erg/(masses[za]+1)
            except KeyError:
                stdout.write(f"WARNING: Capture Collision Missing Mass, {nps} {za}\n")
                dep = prev_erg
        # Alpha production
        elif rxn == 107:
            try:
                dep = prev_erg + qvals[za][rxn]
            except KeyError:
                stdout.write(f"WARNING: Alpha Production Missing Q-value, {nps} {za} {rxn}\n")
                dep = prev_erg
        # Other unhandled reactions
        else:
            stdout.write(f"WARNING: Other Reaction Not Handled, {nps} {za} {prev_erg} {rxn}\n")
            dep = prev_erg-erg
    # --- Photon collisions ---
    elif ipt == 2:
        za = int(za/1000)
        # Coherent/Incoherent scatter: energy difference
        if rxn == -1 or rxn == -2:
            dep = prev_erg-erg
        # Pair production: energy difference minus electron-positron rest mass
        elif rxn == -4:
            dep = prev_erg-1.022
        # Other unhandled reactions
        else:
            stdout.write(f"WARNING: Other Reaction Not Handled, {nps} {za} {prev_erg} {rxn}\n")
            dep = prev_erg-erg
    else:
        dep = 0.0

    return dep

def convert_mt2ntyn(rxn,nps):
    """Converts an MCNP MT reaction number to a PoliMi NTYN code."""
    try:
        return mt2ntyn[rxn]
    except KeyError:
        stdout.write(f"WARNING: MT to NTYN Conversion Failed for rxn={rxn}, nps={nps}\n")
        return rxn # Return original if not found

def all_events_par_and_gen_numbers(event_type, current_particle_num, current_generation,
                                   source_event, bank_event, current_event_data):
    """Updates particle and generation numbers when processing all event types."""
    if event_type == 1000: # Source event
        current_particle_num += 1
        if source_event is not None:
            bank_event = None
            current_generation = 0
        source_event = np.array([
            current_event_data['x'], current_event_data['y'], current_event_data['z'], current_event_data['source_type']
        ])
    elif event_type == 2000: # Bank event
        current_particle_num += 1
        bank_event = np.array([
            current_event_data['x'], current_event_data['y'], current_event_data['z'], current_event_data['bank_type']
        ])
        if source_event is not None and np.sum(bank_event == source_event) != 4:
            current_generation += 1
    return event_type, current_particle_num, current_generation, source_event, bank_event

def surface_crossing_par_number_update(event_type, current_particle_num, current_event_data,
                                       prev_event_data, prev_surface_crossing):
    """Updates particle number based on surface crossing events."""
    if event_type == 3000: # Surface crossing
        if prev_event_data is None or prev_surface_crossing is None:
            prev_surface_crossing = current_event_data
            return current_particle_num + 1, prev_surface_crossing
        
        time_increased = current_event_data['time'] > prev_surface_crossing['time']
        collisions_increased = current_event_data['num_collisions_this_branch'] > prev_surface_crossing['num_collisions_this_branch']
        
        
        if not (time_increased and collisions_increased):
            prev_surface_crossing = current_event_data
            return current_particle_num + 1, prev_surface_crossing
            
    return current_particle_num, prev_surface_crossing

def process_chunk(args):
    """
    Worker function to process a specific chunk of the RecordLog.
    This function is executed by each parallel process.
    """
    ptrac_hdf5_path, start_idx, end_idx, cells, all_events_bool, nEthresh, pEthresh = args
    
    output_lines = []
    
    event_data_lookup = {
        1000: 'Source', 2000: 'Bank', 3000: 'SurfaceCrossing',
        4000: 'Collision', 5000: 'Termination'
    }

    # Each worker opens its own read-only handle to the HDF5 file
    with h5py.File(ptrac_hdf5_path, 'r') as hf:
        record_log_dataset = hf['/ptrack/RecordLog']
        log_chunk = record_log_dataset[start_idx:end_idx]

        # State variables, scoped to this chunk
        last_nps = -1
        current_event_data = None
        
        for rlog_entry in log_chunk:
            nps = rlog_entry['nps']
            event_type = rlog_entry['type']
            event_array_index = rlog_entry['event_array_index']
            
            prev_event_data = current_event_data

            if nps != last_nps:
                current_generation = 0
                current_particle_num = 0
                last_nps = nps
                source_event = None
                bank_event = None
                prev_surface_crossing = None
                last_known_energy = 0.0

            dataset_name = f'/ptrack/{event_data_lookup.get(event_type)}'
            if dataset_name in hf:
                current_event_data = hf[dataset_name][event_array_index]
            else:
                continue # Skip if dataset for event type is missing

            if current_event_data is None:
                continue

            current_particle_energy_after_event = current_event_data['energy']

            if all_events_bool:
                _, current_particle_num, current_generation, source_event, bank_event = \
                    all_events_par_and_gen_numbers(event_type, current_particle_num, current_generation,
                                                   source_event, bank_event, current_event_data)
            else:
                current_particle_num, prev_surface_crossing = \
                    surface_crossing_par_number_update(event_type, current_particle_num, current_event_data,
                                                       prev_event_data, prev_surface_crossing)

            if event_type == 4000: # Collision event
                cell_id = current_event_data['cell_id']
                if (cells is None) or (cell_id in cells):
                    energy_after_collision = current_particle_energy_after_event
                    energy_before_collision = last_known_energy
                    
                    # Simplified energy deposition calculation
                    energy_deposited = energy_before_collision - energy_after_collision

                    if energy_deposited > nEthresh:
                        x, y, z = current_event_data['x'], current_event_data['y'], current_event_data['z']
                        time = current_event_data['time']
                        reaction_code_num = current_event_data['reaction_type']
                        zaid = current_event_data['zaid']
                        particle_type_num = current_event_data['particle_type']
                        weight = current_event_data['weight']
                        num_collisions = current_event_data['num_collisions_this_branch']
                        
                        gen_num = current_generation
                        part_num = current_particle_num
                        
                        reaction_code_num = convert_mt2ntyn(reaction_code_num, nps)
                        
                        output_line = (
                            f"{str(nps):>7s} {str(part_num):>4s} {str(particle_type_num):>2s} "
                            f"{str(reaction_code_num):>4s} {str(zaid):>5s} {str(cell_id):>7s} "
                            f"{energy_deposited:>10.6f} {time:>10.2f} {x:>5.2f} {y:>5.2f} {z:>5.2f} "
                            f"{weight:>11.3E} {str(gen_num):>4s} {str(num_collisions):>4s} {str(0):>4s} "
                            f"{energy_after_collision:>11.3E}\n"
                        )
                        output_lines.append(output_line)

            last_known_energy = current_particle_energy_after_event
            
    return output_lines


def process_ptrac_hdf5_to_polimi_ascii_parallel(ptrac_hdf5_path, polimi_ascii_path,
                                                numCoinc=1, cells=None,
                                                nEthresh=0.0, pEthresh=0.0,
                                                all_events_bool=1,
                                                histories_per_chunk=10000): # Number of HISTORIES per chunk
    """
    Processes MCNP6.3 PTRAC HDF5 output in parallel to a PoliMi ASCII format.
    It parallelizes the processing by distributing independent particle histories (nps)
    across multiple CPU cores.
    """
    print('--- Parallel PTRAC HDF5 to PoliMi Converter ---')
    
    try:
        with h5py.File(ptrac_hdf5_path, 'r') as hf:
            if '/ptrack/RecordLog' not in hf:
                print("Error: '/ptrack/RecordLog' not found.")
                return
            
            print("Scanning RecordLog to identify history boundaries...")
            record_log = hf['/ptrack/RecordLog']
            nps_column = record_log['nps'][:] # Load into memory for faster processing
            
            # Find the indices where the nps value changes
            # This marks the start of a new history's events
            change_indices = np.where(np.diff(nps_column) != 0)[0] + 1
            
            # Create a list of start/end row indices for each history
            history_starts = np.insert(change_indices, 0, 0)
            
            tasks = []
            # *** MODIFIED LOGIC: Chunk by number of histories, not rows ***
            for i in range(0, len(history_starts), histories_per_chunk):
                start_idx = history_starts[i]
                
                # Determine the end index for this chunk
                # It's the start of the next chunk, or the end of the file
                end_idx_pos = i + histories_per_chunk
                if end_idx_pos < len(history_starts):
                    end_idx = history_starts[end_idx_pos]
                else:
                    end_idx = len(nps_column)
                
                args = (ptrac_hdf5_path, start_idx, end_idx, cells, all_events_bool, nEthresh, pEthresh)
                tasks.append(args)

        print(f"Divided work into {len(tasks)} chunks, ensuring complete histories in each.")
        
        with open(polimi_ascii_path, 'w') as f_out:
            # f_out.write(
            #     "History_number Particle_number Particle_Type Reaction_Type Target_nucleus "
            #     "Cell_no._of_collision_event Energy_deposited_in_collision(MeV) "
            #     "Time(shakes) Collision_position_x(cm) Collision_position_y(cm) "
            #     "Collision_position_z(cm) Weight Generation_number "
            #     "Num_collisions_this_branch Code Energy(MeV)\n"
            # )
            
            # Use all available CPU cores
            num_workers = os.cpu_count()
            print(f"Starting parallel processing on {num_workers} cores...")
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
                # Use tqdm for a progress bar
                results = list(tqdm(executor.map(process_chunk, tasks), total=len(tasks), desc="Processing Chunks"))

            print("Writing results to output file...")
            for line_list in results:
                for line in line_list:
                    f_out.write(line)

        print(f"\nSuccessfully processed data to {polimi_ascii_path}")

    except FileNotFoundError:
        print(f"Error: The file '{ptrac_hdf5_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# This guard is crucial for multiprocessing to work correctly
if __name__ == "__main__":
    # Example Usage:
    process_ptrac_hdf5_to_polimi_ascii_parallel(
        ptrac_hdf5_path='timing_and_par_num_ver/U-235/ptrac/two_stilbene/nearly_crit/sur_col/two_stilbene.p.h5',
        polimi_ascii_path='timing_and_par_num_ver/U-235/ptrac/two_stilbene/nearly_crit/sur_col/polimi_collision_data.txt',
        cells=[901,902],
        all_events_bool=0,
	nEthresh=0.03,
        histories_per_chunk=5000 # Adjust this based on memory and number of histories
    )
