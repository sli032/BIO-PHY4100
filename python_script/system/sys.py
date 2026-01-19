import hoomd
from hoomd import md
import json
import copy
import numpy as np
import gsd.hoomd

class System(object):
    def __init__(self, sim, mem, param):
        self.sim = sim
        self.mem = mem
        self.param = copy.copy(param)
        self.system_types = ['M', 'P', 'H']
        self.particlesizes = [self.mem.diameter] * 3
        assert(len(self.system_types) == len(self.particlesizes))
        self.nl = hoomd.md.nlist.Cell(buffer = 0.4, exclusions = ('body', 'constraint'))
        self.bondflip = None


    def dump_bonds(self):
        if self.bondflip:
            a = self.bondflip.get_bond_neigh_bonds()
            b = self.bondflip.get_bond_dihedrals()
            c = self.bondflip.get_bond_faces()
        else:
            a = self.mem.bond_neigh_bonds
            b = self.mem.bond_dihedrals
            c = self.mem.bond_faces
        f = open("bond_neigh_bonds.json", "w")
        json.dump(a, f)
        f.close()
        f = open("bond_dihedrals.json", "w")
        json.dump(b, f)
        f.close()
        f = open("bond_faces.json", "w")
        json.dump(c, f)
        f.close()

    
    def setup_potential(self):
        self.lj = hoomd.md.pair.LJ(self.nl, default_r_cut = self.param['lj_cut'])
        self.lj.mode = 'shift'

        self.bond_harmonic = hoomd.md.bond.Harmonic()
        
        self.dihedral_force = hoomd.md.dihedral.Harmonic()
        if self.param['bondflip']:
            self.bondflip = hoomd.bondflip_plugin.update.bondflip(hoomd.trigger.Periodic(self.param['bondflip_period']),
            self.mem.bond_neigh_bonds, self.mem.bond_dihedrals, 
            self.mem.bond_faces,
            bond_energy_name = 'harmonic',
            k = self.param['mbondk'], 
            r0 = self.param['mbond_r0'], 
            delta = 0.,
            dihedral = True, 
            angle = False,
            kappa = self.param['mdihedk'], 
            phi_0 = np.pi-self.param['mtheta0'], 
            T = self.param['bondflip_kT'])

    def setup_potential_coef(self):
        # lj potential
        for i in range(len(self.system_types)):
            for j in range(len(self.system_types)):
                if i <= j:
                    type1, size1 = self.system_types[i], self.particlesizes[i]
                    type2, size2 = self.system_types[j], self.particlesizes[j]
                    dis = 0.5*(size1+size2)
                    lj = self.param.get('lj_'+type1+type2, 0.) + self.param.get('lj_'+type2+type1, 0.)
                    if type1 == type2:
                        lj = 0.5 * lj
                    print(type1, type2, dis, 'interact', lj)
                    if lj > 0.:
                        print(type1, type2, dis, 'interact', lj)
                        self.lj.params[(type1, type2)] = dict(sigma = dis*2**(-1/6), epsilon = lj)
                        self.lj.r_cut[(type1, type2)] = self.param['lj_cut']
                    else:
                        print(type1, type2, dis, 'exclude', lj)
                        self.lj.params[(type1, type2)] = dict(sigma = dis*2**(-1/6), epsilon = self.param['lj_repel'])
                        self.lj.r_cut[(type1, type2)] = dis
                    
        # bond
        self.bond_harmonic.params['mbond'] = dict(k = self.param['mbondk'], r0 = self.param['mbond_r0'])
        
        # dihedral
        self.dihedral_force.params['mdihedral'] = dict(k = 2 * self.param['mdihedk'], d = -1, n = 1, phi0 = np.pi-self.param['mtheta0'])
    
    def setup_operation(self):
        gsd_writer = hoomd.write.GSD(filename="md.gsd", filter=hoomd.filter.All(), trigger=hoomd.trigger.Periodic(self.param['record_period']), mode='ab', dynamic=['property', 'topology'])
        self.sim.operations.writers.append(gsd_writer)

        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All())
        self.sim.operations.computes.append(thermodynamic_properties)

        logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        logger.add(thermodynamic_properties, quantities=['potential_energy'])
        logger.add(self.sim, quantities=['timestep', 'walltime'])
        logger.add(self.lj, quantities=['energy'])
        logger.add(self.bond_harmonic, quantities=['energy'])
        logger.add(self.dihedral_force, quantities=['energy'])
        if self.bondflip:
            logger.add(self.bondflip, quantities=['num_accepted_bonds'])
            logger.add(self.bondflip, quantities=['accepted_bonds'])
        file = open('log.txt', mode='a', newline='\n')
        table_file = hoomd.write.Table(output=file,
            trigger=hoomd.trigger.Periodic(self.param['record_period']),
            logger=logger)
        self.sim.operations.writers.append(table_file)

        if self.bondflip:
            self.sim.operations.updaters.append(self.bondflip)

    def setup_integrator(self):
        langevin = hoomd.md.methods.Langevin(filter = hoomd.filter.All(), kT = self.param['kT'])

        self.sim.operations.integrator = hoomd.md.Integrator(dt = self.param['dt'], methods = [langevin], forces = [self.lj, self.bond_harmonic, self.dihedral_force])
    
    def gsd_writer(self, snap = None):
        if snap is None:
            snap = self.sim.state.get_snapshot()
        try:
            frame = gsd.hoomd.Frame()
            frame.configuration.box = snap.configuration.box
            frame.particles.N = snap.particles.N
            frame.particles.types = snap.particles.types
            frame.particles.typeid = [
                1 if len(snap.bonds.group[np.where(snap.bonds.group == i)]) == 5 else
                2 if len(snap.bonds.group[np.where(snap.bonds.group == i)]) == 6 else
                0
                for i in range(snap.particles.N)]
            frame.particles.position = snap.particles.position
            frame.particles.diameter = snap.particles.diameter
            frame.particles.velocity=snap.particles.velocity
            frame.bonds.N = snap.bonds.N
            frame.bonds.types = snap.bonds.types
            frame.bonds.typeid = snap.bonds.typeid
            frame.bonds.group = snap.bonds.group
            frame.dihedrals.N = snap.dihedrals.N
            frame.dihedrals.types = snap.dihedrals.types
            frame.dihedrals.typeid = snap.dihedrals.typeid
            frame.dihedrals.group = snap.dihedrals.group
            with gsd.hoomd.open(name='color.gsd', mode='a') as f:
                f.append(frame)
        except Exception as e:
            print(f"Error writing to GSD file: {e}", 'red')