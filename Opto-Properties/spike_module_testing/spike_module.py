from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from neuron import h
import time

pc = h.ParallelContext()
MPI_RANK = int(pc.id())
N_HOSTS = int(pc.nhost())

class Spiker(SimulatorMod):
    def __init__(self, node_set=None):
        """Initialize the Spiker module."""
        self._spike_records = {}

    def _set_spike_detector(self, sim):
        for gid, cell in sim.net.get_local_cells().items():
            tvec = sim._spikes[gid]
            self._spike_records[gid] = tvec

    def initialize(self, sim):
        """Called once at the beginning of the simulation."""
        self._set_spike_detector(sim)

    def block(self, sim, block_interval):
        """Called every block to collect and print spike data."""
        
        # Use py_gather to collect spike data from all ranks on rank 0
        gathered_data = pc.py_gather(self._spike_records, 0)
        pc.barrier()

        if MPI_RANK == 0:
            # Combine all rank_data into one dictionary
            combined_spike_data = {}
            for rank_data in gathered_data:
                combined_spike_data.update(rank_data)
            
            # Sort the dictionary by gid not needed but makes the printout of spikes look nice 
            sorted_spike_data = dict(sorted(combined_spike_data.items()))

            # Print the sorted spike data
            for gid, spikes in sorted_spike_data.items():
                print(f"spikes for gid {gid}: {spikes.to_python()} (len={len(spikes)})")
        
        # Set the spike detector for the next block
        self._set_spike_detector(sim)

    def finalize(self, sim):
        """Called at the end of the simulation."""
        # Any finalization code you need
        pass
