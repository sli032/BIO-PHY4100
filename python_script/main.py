import sys
import hoomd
import hoomd.md
import hoomd.bondflip_plugin
import system
import numpy as np
import os
import json, glob

    

def init(file):
    cpu = hoomd.device.CPU()

    snap = hoomd.Snapshot()
    snap.configuration.box = [param['boxL'], param['boxL'], param['boxL'], 0, 0, 0]
    snap.particles.types = ['M', 'P', 'H']
    snap.bonds.types = ['mbond']
    snap.dihedrals.types = ['mdihedral']
    typedic = {'M':0, 'P':1, 'H':2, 'mbond':0, 'mangle':0, 'mdihedral':0}

    mem = system.membrane.mem(param)
    mem.frame_in()
    mem.update_mem_info()
    mem.update_bond_flip_info()
    system.snap.update_snap(snap, mem, typedic, bond = True, angle = False, dihedral = True)

    sim = hoomd.Simulation(device = cpu, seed = param['seed'])
    sim.create_state_from_snapshot(snap)
    memsys = system.sys.System(sim, mem, param)
    print('sim.particles.types:', memsys.sim.state.particle_types)

    memsys.setup_potential()
    memsys.setup_potential_coef()
    memsys.setup_operation()
    memsys.setup_integrator()

    return memsys

def init_from_file(file):
    cpu = hoomd.device.CPU()

    print("continue running")
    
    mem = system.membrane.mem(param)
    if param['bondflip']:
        mem.update_bond_from_file()
    
    sim = hoomd.Simulation(device = cpu, seed = param['seed'])
    sim.create_state_from_gsd(file, frame = -1)
    memsys = system.sys.System(sim, mem, param)

    memsys.setup_potential()
    memsys.setup_potential_coef()
    memsys.setup_operation()
    memsys.setup_integrator()

    return memsys


def run(sys): 
    time = 0    
    while time < param['runtime']:
        sys.sim.run(param['record_period'], write_at_start= True)
        sys.gsd_writer()
        time += param['record_period']
    sys.dump_bonds()

def anneal_run(sys): 
    while param['kT'] >= 0.:            
        print('running with temp ', param['kT'], 'bondflip_kT ', param['bondflip_kT'])
        time = 0
        while time < param['runtime']:
            sys.sim.run(param['record_period'], write_at_start= True)
            sys.gsd_writer()
            time += param['record_period']
        sys.dump_bonds()
        param['kT'] = round(param['kT'] - 0.02, 2)
        param['bondflip_kT'] = round(param['bondflip_kT'] - 0.02, 2)
        sys = init_from_file(param['gsdfile']+'.gsd')



param=sys.argv[1]
print(param)
with open(param) as f: param=json.load(f)

if glob.glob(param['gsdfile']+'.gsd'):
    memsys = init_from_file(param['gsdfile']+'.gsd')
else:
    memsys = init(param['gsdfile']+'.gsd')
if param['annealrun']:
    anneal_run(memsys)
else:
    run(memsys)