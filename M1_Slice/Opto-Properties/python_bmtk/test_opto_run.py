from bmtk.simulator import bionet
from opto_stim import optoStim
import random
import numpy as np


random.seed(10)
np.random.seed(10)
conf = bionet.Config.from_json('config.json', validate=True)
conf.build_env()

stim_params = {
    "nPulses": 1,
    "Dt_on": 2,
    "Dt_off": 95,
    "gmax": 0.04,
    "tauChR2": 0.4,
    "Gd1": 0.25,
    "Gd2": 0.5,
    "light_intensity": 1,
}

opto_module = optoStim(
    node_set="PN_cells",
    section_name="soma",
    stim_params=stim_params,
    pulse_file="input_test.csv",
    mat_file_path='LightData_blue.mat',
    probe_file_path='probe_positions_full.csv',
    uptake_prob=1)

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)
sim.add_mod(opto_module)  
sim.run()
bionet.nrn.quit_execution()
