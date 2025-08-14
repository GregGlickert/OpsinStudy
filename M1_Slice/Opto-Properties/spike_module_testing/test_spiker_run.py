from bmtk.simulator import bionet
from opto_stim import optoStim
from spike_module import Spiker


conf = bionet.Config.from_json('config.json', validate=True)
conf.build_env()


spikeModule = Spiker()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)
sim.add_mod(spikeModule)  
sim.run()
bionet.nrn.quit_execution()
