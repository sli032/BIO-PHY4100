# Copyright (c) 2009-2022 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""BondFlip Updater."""

# Import the C++ module.
from hoomd.bondflip_plugin import _bondflip_plugin
# import copy
from hoomd.logging import log
from hoomd.data.parameterdicts import ParameterDict
# from hoomd.logging import Loggable
from hoomd.operation import Updater


# Import the hoomd Python package.
import hoomd 


class bondflip(Updater):
    # Initialize the bondflip plugin
    #
    # \param bond_neigh_bonds : Dictionary of bond index and its four neighbor edges' indices.
    #        bond_dihedral: Dictionary of bond index and its dihedral index.
    #        bond_faces: Dictionary of bond index and its triangle indices.
    #        group : Group used to calculate the temperature.
    #        period: Time steps between two bondflip attempts.
    # def __init__(self, bond_neigh_bonds, bond_dihedral, bond_faces):
    def __init__(self,  
                trigger, 
                bond_neigh_bonds, 
                bond_dihedral, 
                bond_faces,
                bond_energy_name = "harmonic",
                k = 1.0,
                r0 = 1.0,
                sigma = 0.0,
                epsilon = 0.0,
                delta = 0.0,
                dihedral = False,
                angle = False,
                kappa = 1.0,
                phi_0 = 0.0,
                T = 1.0,
                Na = 4,
                Nb = 7):
        
        print("initializing bondflip in python")

        self.bond_neigh_bonds = bond_neigh_bonds
        self.bond_dihedral = bond_dihedral
        self.bond_faces = bond_faces

        params = ParameterDict(
            bond_energy_name = "harmonic",
            k = 1.0,
            r0 = 1.0,
            sigma = 0.0,
            epsilon = 0.0,
            delta = 0.0,
            dihedral = False,
            angle = False,
            kappa = 1.0,
            phi_0 = 0.0,
            T = 1.0,
            Na = 4,
            Nb = 7)
        params.update(dict(
            bond_energy_name = bond_energy_name,
            k = k,
            r0 = r0,
            sigma = sigma,
            epsilon = epsilon,
            delta = delta,
            dihedral = dihedral,
            angle = angle,
            kappa = kappa,
            phi_0 = phi_0,
            T = T,
            Na = Na,
            Nb = Nb))
        self._param_dict.update(params)

        print(self._param_dict)
        # This updater has to be applied every timestep
        super().__init__(trigger)

    def _attach(self):
        # if isinstance(self._simulation.device, hoomd.device.CPU):
        self._cpp_obj = _bondflip_plugin.BondFlipUpdater(
                            self._simulation.state._cpp_sys_def, 
                            self.trigger,
                            self.bond_neigh_bonds, 
                            self.bond_dihedral, 
                            self.bond_faces)

        self._cpp_obj.set_params(self.bond_energy_name, self.k, self.r0, self.sigma, self.epsilon, self.delta, self.dihedral, self.angle, self.kappa, self.phi_0, self.T, self.Na, self.Nb)
        super()._attach()

    @log(category="scalar", requires_run=True)
    def num_accepted_bonds(self):
        return self._cpp_obj.get_num_accepted_bonds()
    
    @log(category="string", requires_run=True)
    def accepted_bonds(self):
        return self._cpp_obj.get_accepted_bonds()

    def get_bond_neigh_bonds(self):
        return self._cpp_obj.get_bond_neigh_bonds()
    
    def get_bond_dihedrals(self):
        return self._cpp_obj.get_bond_dihedrals()
    
    def get_bond_faces(self):
        return self._cpp_obj.get_bond_faces()
