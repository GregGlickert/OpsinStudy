from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from bmtk.simulator.bionet.io_tools import io
from neuron import h
import pandas as pd 
import numpy as np
import scipy.io
import h5py
import os

# Parallel Context to manage distributed execution across MPI ranks
pc = h.ParallelContext()
MPI_RANK = int(pc.id())  # Unique rank ID for this process
N_HOSTS = int(pc.nhost())  # Total number of processes (hosts)

class optoStim(SimulatorMod):
    def __init__(self, node_set, section_name, section_index=0, section_dist=0.5, stim_params=None, pulse_file=None,output_file="opto_stim_currents.h5", tmp_dir="./tmp", mat_file_path=None, probe_file_path=None, uptake_prob=1.):
        """
        Initialize the optoStim module for applying optogenetic stimulation.
        
        :param node_set: Name or ID of the node set to apply the stimulation.
        :param section_name: The section of the cell to apply stimulation (e.g., 'soma').
        :param section_index: Index of the section in the structure (default is 0).
        :param section_dist: The location along the section (0.0 to 1.0, default is 0.5).
        :param stim_params: Dictionary containing the stimulation parameters (e.g., light intensity, duration).
        :param pulse_file: Path to a file containing pulse times (one pulse per line).
        :param output_file: HDF5 file to save the output currents (default is 'opto_stim_currents.h5').
        :param tmp_dir: Directory to store temporary files during simulation.
        """
        # Initialize instance variables based on input parameters
        self._node_set = node_set
        self._section_name = section_name
        self._section_index = section_index
        self._section_dist = section_dist
        self._stim_params = stim_params or {}
        self._pulse_file = pulse_file
        self._pulse_times = None
        self._output_file = output_file if os.path.isabs(output_file) else os.path.join(tmp_dir, output_file)
        self._tmp_dir = tmp_dir
        self._tmp_file = self._get_tmp_fname(MPI_RANK)
        self._tmp_handle = None
        self._block_data = None
        self._gids_on_rank = None
        self.mat_file_path = mat_file_path
        self.probe_file_path = probe_file_path
        self.uptake_prob = uptake_prob

        # If a pulse file is provided, load the pulse times
        if pulse_file:
            self._load_pulse_file(pulse_file)

        # Initialize lists for storing stimulation and connection data
        self._stims = []
        self._netcons = []
        self._vecstims = []
        self.currents = {}

        # If the output file already exists, remove it to avoid overwriting errors
        if MPI_RANK == 0:
            if os.path.isfile(self._output_file):
                os.remove(self._output_file)
        pc.barrier() # wait up guys dont leave without me 


    def _get_tmp_fname(self, rank):
        """Generate a unique temporary filename for each rank."""
        return os.path.join(self._tmp_dir, f'tmp_{rank}_{os.path.basename(self._output_file)}')


    def _create_tmp_file(self, nsteps, nsites, nprobes):
        """Create a temporary HDF5 file for storing current data during the simulation."""
        os.makedirs(self._tmp_dir, exist_ok=True)
        self._tmp_handle = h5py.File(self._tmp_file, 'w')
        
        # Create a dataset to store the stimulation data for each time step, site, and probe
        self._tmp_handle.create_dataset(
            '/stim/data', (nsteps, nsites, nprobes), dtype=np.float64, chunks=True, maxshape=(None, nsites, nprobes)
        )
        
        # Create a dataset to store the GIDs of the sites (cells)
        self._tmp_handle.create_dataset('/stim/gids', (nsites,), dtype=np.int32)  
        self._tmp_handle['/stim/gids'][:] = self._gids_with_virus  # Save GIDs for each rank
        self._tmp_handle.flush()


    def _write_block_data(self, block_interval):
        """Write current data to the temporary file at each block interval."""
        start_step = block_interval[0]
        end_step = block_interval[1]
        
        # Store the data for the specified block range in the temporary file
        self._tmp_handle['/stim/data'][start_step:end_step, :, :] = self._block_data[start_step:end_step, :, :]
        self._block_data[start_step:end_step, :, :] = 0  # Reset the block data after writing
        self._tmp_handle.flush()


    def _combine_files(self, nsteps, total_nsites,nprobes):
        """Combine all temporary files from each MPI rank into the final output file."""
        if MPI_RANK == 0:
            with h5py.File(self._output_file, 'w') as out_file:
                # Create a combined dataset in the output file to store data from all ranks
                out_file.create_dataset(
                    '/stim/data', (nsteps, total_nsites, nprobes), dtype=np.float64, chunks=True
                )
                out_file.create_dataset('/stim/gids', (total_nsites,), dtype=np.int32)  # Combined GIDs

                # Initialize a list to track the global index for each rank's data
                current_idx = 0

                # Iterate over each rank to collect and combine data
                for rank in range(N_HOSTS):
                    tmp_file = self._get_tmp_fname(rank)
                    with h5py.File(tmp_file, 'r') as tmp_handle:
                        # Get the number of sites (cells) on this rank
                        nsites_on_rank = len(tmp_handle['/stim/gids'][:])
                        # Copy data for the current rank into the final output file
                        out_file['/stim/data'][:, current_idx:current_idx + nsites_on_rank, :] = tmp_handle['/stim/data'][:]
                        out_file['/stim/gids'][current_idx:current_idx + nsites_on_rank] = tmp_handle['/stim/gids'][:]

                        # Update the global index to point to the next segment
                        current_idx += nsites_on_rank
        pc.barrier()



    def initialize(self, sim):
        """Called once at the start of the simulation to set up the stimulation."""
        io.log_info("Setting up optogenetic stimulation")
        io.log_info(f"Probability of infection is {self.uptake_prob}")
        # Get the GIDs for the selected node set
        select_gids = list(sim.net.get_node_set(self._node_set).gids())
        
        # Determine the local GIDs assigned to this rank
        self._gids_on_rank = list(set(select_gids) & set(sim.local_gids))
        
        # Decide which cells receive the virus based on uptake probability
        self._gids_with_virus = []
        for gid in self._gids_on_rank:
            if np.random.rand(1) < self.uptake_prob:
                self._gids_with_virus.append(gid)
                
        nsteps = int(sim.n_steps)  # Number of simulation time steps
        
        # get the gids that will be infected
        
        nsites = len(self._gids_with_virus)  # Number of sites (cells) for this rank
        
        # Read the probe file to get the probe positions
        probe_df = pd.read_csv(self.probe_file_path)
        nprobes = len(probe_df)  # Number of probes
        io.log_info(f"Number of light sources is {nprobes}")
        
        # Create temporary file to store current data
        self._create_tmp_file(nsteps, nsites, nprobes)
        
        # Initialize block data array to store current values
        self._block_data = np.zeros((nsteps, nsites, nprobes))

        # Loop through the GIDs to apply the stimulation
        for gid in self._gids_with_virus:
            self.currents[gid] = []
            for index, row in probe_df.iterrows():
                probe_pos = np.array([row['X'], row['Y'], row['Z']])
                probe_number = row['probe_number']
                #print(f"probe location is {probe_pos}")
                
                # Get the cell corresponding to the current GID
                cell = sim.net.get_cell_gid(gid)
                
                # Get the section of the cell where stimulation will be applied
                hobj_sec = getattr(cell.hobj, self._section_name)[self._section_index](self._section_dist)
                
                # Create the stimulation for this section and cell
                self._create_stim(hobj_sec, cell, probe_pos, gid, probe_number)


    def step(self, sim, tstep):
        """Called at each simulation time step to update the stimulation currents."""
        for gid, vec in self.currents.items():
            # Find the local index of this GID on the current rank
            local_idx = self._gids_with_virus.index(gid)

            # Store current value for each probe at the current time step
            for j, probe in enumerate(vec):
                #print(f"Storing current for gid {gid}, probe {j}, tstep {tstep}")
                self._block_data[tstep, local_idx, j] = probe.x[-1]  # Store at the correct index


    def block(self, sim, block_interval):
        """Called every block."""
        self._write_block_data(block_interval)


    def finalize(self, sim):
        """Called at the end of the simulation."""
        # Close the temporary file handle
        self._tmp_handle.close()

        # Sync across all ranks
        pc.barrier()

        # Number of time steps (simulation steps)
        nsteps = sim.n_steps

        # Number of sites (cells) assigned to the current rank
        local_nsites = len(self._gids_with_virus)  # Local number of infected cells on this rank
        total_nsites = pc.allreduce(local_nsites, 1) # 1 will sum all local_nsites and give each rank the total number of sites

        # Read the probe file to determine the number of probes
        probe_df = pd.read_csv(self.probe_file_path)
        nprobes = len(probe_df)

        # Combine temporary files into the final output file, considering sites and probes
        self._combine_files(nsteps, total_nsites, nprobes)

        # Clean up the temporary files on all ranks
        if MPI_RANK == 0:
            for rank in range(N_HOSTS):
                tmp_file = self._get_tmp_fname(rank)
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)


    ## CODE FOR OPTO SETUP ##
    def _load_pulse_file(self, pulse_file):
        """Load pulse times from a text file (one pulse per line)."""
        try:
            self._pulse_times_df = pd.read_csv(pulse_file)
        except Exception as e:
            raise ValueError(f"Error loading pulse file '{pulse_file}': {e}")
    
    
    def _get_light_scattering(self,mat):
        """Reads in the light scattering data from the monte carlo simulation

        Args:
            mat (mat): datafile from simulation

        """
        mat = scipy.io.loadmat(mat)
        # Extract the data arrays
        r = mat['r'].flatten()  # Radial distances (r)
        depth = mat['depth'].flatten()  # Depth (depth)
        plot_mat = mat['plot_mat']  # The matrix to plot
        
        r = r * 1000  # Convert radial distance to micrometers (µm)
        depth = depth * 1000  # Convert depth to micrometers (µm)

        # Create mirrored radial distance (r_mirrored)
        r_mirrored = np.concatenate([-r[::-1], r])  # Create a new range from -r to +r

        # Mirror the plot_mat matrix: reverse the columns of plot_mat and concatenate
        plot_mat_mirrored = np.concatenate([plot_mat[:, ::-1], plot_mat], axis=1)

        # Replace -inf values with np.nan
        plot_mat_mirrored[plot_mat_mirrored == -np.inf] = np.nan

        # Calculate the min and max ignoring NaN
        min_val = np.nanmin(plot_mat_mirrored)
        max_val = np.nanmax(plot_mat_mirrored)

        plot_mat_mirrored_normalized = (plot_mat_mirrored - min_val) / (max_val - min_val)

        # Create a meshgrid for the mirrored r and depth
        R, Depth = np.meshgrid(r_mirrored, depth)

        return R,Depth,plot_mat_mirrored_normalized
    
    
    def _get_intensity(self,cell,R,Depth,plot_mat,probe_position=[0,0,0]):
        """gets the light intensity for every node.

        Args:
            cell 
            R (2d np array): radidal distance from monte carlo sim
            Depth (2d np array): depth from monte carlo sim
            plot_mat (2d np array): light intensity with R and Depth from monte carlo 
            probe_position (tuple): where light comes out from. Defaults to [0,0,0]
        """
        cell_location = np.array(cell.soma_position)

        distance_x = cell_location[0] - probe_position[0]
        distance_y = cell_location[1] - probe_position[1]
        distance_z = cell_location[2] - probe_position[2]

        radial_distance = np.sqrt(distance_x**2 + distance_y**2)

        r_indices = np.array([np.argmin(np.abs(R[0, :] - r)) for r in radial_distance])
        depth_indices = np.array([np.argmin(np.abs(Depth[:, 0] - z)) for z in distance_z])


        node_intensities = plot_mat[depth_indices, r_indices]
        if np.any(np.isnan(node_intensities)):
            print(f"Warning: NaN values detected in node intensities at indices {np.isnan(node_intensities)}")

        node_intensities = np.nan_to_num(node_intensities, nan=0)
        

        return node_intensities


    def _create_stim(self, hobj_sec, cell,probe_pos,gid,probe_number):
        """Create the optogenetic stimulation for the specified section."""
        stim = h.ChR2_william_event(hobj_sec)
        
        max_intensity = self._stim_params.get("light_intensity")
        rad_dist,depth_dist,intensity_dist = self._get_light_scattering(self.mat_file_path)
        norm_intensity = self._get_intensity(probe_position=probe_pos,cell=cell,R=rad_dist,Depth=depth_dist,plot_mat=intensity_dist)
        light_intensity = max_intensity*norm_intensity 
        #print(f"Light intensity is {light_intensity} for gid {gid}")
        stim.light_intensity = light_intensity
        
        stim.nPulses = self._stim_params.get("nPulses")
        stim.Dt_on = self._stim_params.get("Dt_on")
        stim.Dt_off = self._stim_params.get("Dt_off")
        stim.gmax = self._stim_params.get("gmax")
        stim.tauChR2 = self._stim_params.get("tauChR2")
        stim.Gd1 = self._stim_params.get("Gd1")
        stim.Gd2 = self._stim_params.get("Gd2")
        
        self._stims.append(stim)

        i = h.Vector().record(stim._ref_i)
        self.currents[gid].append(i)

        pulse_times = self._pulse_times_df[self._pulse_times_df["probe_number"] == probe_number]["timestamps"].values
        #print(f"pulse times are {pulse_times}")
        pulse_vec = h.Vector(pulse_times)
        vec_stim = h.VecStim()
        vec_stim.play(pulse_vec)
        self._vecstims.append(vec_stim)

        netcon = h.NetCon(vec_stim, stim)
        netcon.weight[0] = 1
        netcon.delay = 0
        self._netcons.append(netcon)


    