import sys
import os
import warnings
import synapses
import pathlib
from bmtk.simulator import bionet
from bmtk.simulator.bionet.pyfunction_cache import add_weight_function
from opto_stim import optoStim
from neuron import h
import time
import shutil

CONFIG = 'config.json'
USE_CORENEURON = False

import os
import json

def get_synaptic_params(path_to_syn_folder):
    """Gets values of all json files and puts them into one file"""
    combined_data = []
    # List all files in the input folder
    for filename in os.listdir(path_to_syn_folder):
        if filename.endswith('.json'):
            file_path = os.path.join(path_to_syn_folder, filename)
            
            # Open and read the JSON file
            with open(file_path, 'r') as file:
                data = json.load(file)
                # Append the filename and its data
                combined_data.append({
                    'filename': filename,
                    'data': data
                })
    return combined_data

def save_synaptic_params(data,path_to_output_dir):
    """Saves combined json data into one file"""
    with open(path_to_output_dir, 'w') as output_file:
        json.dump(data, output_file, indent=4)

def run(config_file=CONFIG, use_coreneuron=USE_CORENEURON):

    warnings.simplefilter(action='ignore', category=FutureWarning)
    pc = h.ParallelContext()
    with open(config_file, 'r') as json_file:
        conf_dict = json.load(json_file)
        if os.environ.get("OUTPUT_DIR"):
            output_dir = os.path.abspath(os.environ.get('OUTPUT_DIR'))
            pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
            conf_dict['manifest']['$OUTPUT_DIR'] = output_dir
            synaptic_report_dir = output_dir + "/synaptic_report.json"
        # Handle COMPONENT_PATH
        if os.environ.get("COMPONENT_PATH"):
            com_path =  os.path.abspath(os.environ.get('COMPONENT_PATH'))
            network_config = conf_dict['network']  # Safely fetch key
            if pc.id() == 0: # only let one node read and edit the json
                with open(network_config, 'r') as net_file:
                    network_dict = json.load(net_file)
                    network_dict['manifest']['$COMPONENTS_DIR'] = com_path

                with open(network_config, 'w') as f:
                    json.dump(network_dict,f,
                        indent=4,  # Ensures 4-space indentation
                        ensure_ascii=False,  # Ensures Unicode characters are not escaped
                        sort_keys=False  # Keeps the key order as written
                    )              
            
                syn_data = get_synaptic_params(os.path.join(com_path,'synaptic_models/synapses_STP'))
        else:
            syn_data = get_synaptic_params('components/synaptic_models/synapses_STP')
            
    pc.barrier() 
    # register synaptic weight function
    synapses.load(randseed=1111)
    add_weight_function(synapses.lognormal_weight, name='lognormal_weight')

    if use_coreneuron:
        import corebmtk
        conf = corebmtk.Config.from_json(conf_dict, validate=True)
    else:
        conf = bionet.Config.from_json(conf_dict, validate=True)

    conf.build_env()
    graph = bionet.BioNetwork.from_config(conf)

    if use_coreneuron:
        sim = corebmtk.CoreBioSimulator.from_config(
            conf, network=graph, gpu=False)
    else:
        sim = bionet.BioSimulator.from_config(conf, network=graph)

    '''
    # This calls insert_mechs() on each cell to use its gid as a seed
    # to the random number generator, so that each cell gets a different
    # random seed for the point-conductance noise
    cells = graph.get_local_cells()
    for cell in cells:
        cells[cell].hobj.insert_mechs(cells[cell].gid)
    '''

    # clear ecp temporary directory to avoid errors
    if pc.id() == 0:
        try:
            ecp_tmp = conf['reports']['ecp']['tmp_dir']
        except:
            pass
        else:
            if os.path.isdir(ecp_tmp):
                for f in os.listdir(ecp_tmp):
                    if f.endswith(".h5"):
                        try:
                            os.remove(os.path.join(ecp_tmp, f))
                        except Exception as e:
                            print(f'Failed to delete {f}. {e}')
    pc.barrier()

    with open(config_file) as f:
        config = json.load(f)
        stim_params = config['input']['opto_input']['stim_params']  
        pulse_file=config['input']['opto_input']['pulse_file']
        probe_file_path=config['input']['opto_input']['probe_file_path']
        node_set = config['input']['opto_input']['node_set']
        section_name = config['input']['opto_input']['section_name']
        intensity_file_path = config['input']['opto_input']['intensity_file_path']
        uptake_prob = config['input']['opto_input']['uptake_prob']
        
    opto_module = optoStim(
        node_set=node_set,
        section_name=section_name,
        stim_params=stim_params,
        tmp_dir=output_dir,
        intensity_file_path=intensity_file_path,
        pulse_file=pulse_file,
        probe_file_path=probe_file_path,
        uptake_prob=uptake_prob)

    sim.add_mod(opto_module)  
    sim.run()
    # must be ran after sim.run since that creates dir
    if pc.id() == 0:
        save_synaptic_params(syn_data,synaptic_report_dir)
        path = os.path.split(output_dir)[0]
        if os.environ.get("COMPONENT_PATH"):
            pass
            #shutil.move(os.environ.get('COMPONENT_PATH'),path)
            # could also just delete it with shutil.rmtree 

    bionet.nrn.quit_execution()


if __name__ == '__main__':
    for i, s in enumerate(sys.argv):
        if s in __file__:
            break

    if i < len(sys.argv) - 1:
        argv = sys.argv[i + 1:]
        for i in range(1, len(argv)):
            try:
                argv[i] = eval(argv[i])
            except:
                argv[i] = argv[i]
        run(*argv)
    else:
        run()
