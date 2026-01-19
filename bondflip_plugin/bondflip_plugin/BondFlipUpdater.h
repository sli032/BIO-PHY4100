// Copyright (c) 2009-2022 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

// **********************
// This is a simple example code written for no function purpose other then to demonstrate the steps
// needed to write a c++ source code plugin for HOOMD-Blue. This example includes an example Updater
// class, but it can just as easily be replaced with a ForceCompute, Integrator, or any other C++
// code at all.

// inclusion guard
#ifndef _BONDFLIP_UPDATER_H_
#define _BONDFLIP_UPDATER_H_

// #include "hoomd/ComputeThermo.h"
// #include "hoomd/Logger.h"
#include "hoomd/ParticleGroup.h"
#include <map>
#include <random>

/*! \file BondFlipUpdater.h
    \brief Declaration of BondFlipUpdater
*/

#include <hoomd/Updater.h>

// pybind11 is used to create the python bindings to the C++ object,
// but not if we are compiling GPU kernels
#ifndef __HIPCC__
#include <pybind11/pybind11.h>
#endif

namespace hoomd
    {
// (if you really don't want to include the whole hoomd.h, you can include individual files IF AND
// ONLY IF hoomd_config.h is included first) For example: #include <hoomd/Updater.h>

// Second, we need to declare the class. One could just as easily use any class in HOOMD as a
// template here, there are no restrictions on what a template can do

//! A nonsense particle updater written to demonstrate how to write a plugin
/*! This updater simply sets all of the particle's velocities to 0 when update() is called.
 */
class BondFlipUpdater : public Updater
    {
    public:
        //! Constructor
        BondFlipUpdater(std::shared_ptr<SystemDefinition> sysdef, 
                        std::shared_ptr<Trigger> trigger,
                        // std::shared_ptr<ComputeThermo> thermo,
                        // std::shared_ptr<Logger> logger,
                        std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
                        std::map<unsigned int, unsigned int> bond_dihedral,
                        std::map<unsigned int, std::vector<unsigned int>> bond_faces);
        void set_params(std::string bond_energy_name,
                        Scalar k,
                        Scalar r0,
                        Scalar sigma,
                        Scalar epsilon,
                        Scalar delta,
                        bool dihedral,
                        bool angle,
                        Scalar kappa,
                        Scalar phi_0,
                        Scalar T,
                        int Na,
                        int Nb);

        //! Take one timestep forward
        virtual void update(uint64_t timestep);
        // virtual std::vector< std::string > getProvidedLogQuantities(void);
        // virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);
        Scalar get_num_accepted_bonds();
        std::string get_accepted_bonds();
        std::map<unsigned int, std::vector<unsigned int> > get_bond_neigh_bonds();
        std::map<unsigned int, unsigned int> get_bond_dihedrals();
        std::map<unsigned int, std::vector<unsigned int> > get_bond_faces();
    
    protected:
        std::mt19937 m_rnd;
        unsigned int m_seed;                               
        std::string m_bond_energy_name;
        Scalar m_k = 1.0;
        Scalar m_r0 = 1.0;
        Scalar m_sigma = 0.0;
        Scalar m_epsilon = 0.0;
        Scalar m_delta = 0.0;
        Scalar m_T = 1.0;
        int m_Na = 4;
        int m_Nb = 7;
        bool m_dihedral = false;
        bool m_angle = false;
        Scalar m_kappa = 1.0;
        Scalar m_phi_0 = 3.1415926;
        int m_accepted = 0;

        bool is_connected(unsigned int, unsigned int);
        void find_opposite_vertex(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int &, unsigned int &);
        void flip_bond(unsigned int, Scalar);
        void printcheck(unsigned int);
        
        Scalar distance(BoxDim, ArrayHandle<Scalar4>, unsigned int, unsigned int);
        bool face_are_neighbor(unsigned int f0, unsigned int f1);
        void update_neigh_bonds(unsigned int, unsigned int, std::vector<unsigned int>, std::map<unsigned int, unsigned int>);
        bool get_bond_energy(BoxDim box, ArrayHandle<Scalar4> h_pos, unsigned int idx_a, unsigned int idx_b, Scalar& bond_eng);
        Scalar get_dihedral_energy(BoxDim, ArrayHandle<Scalar4>, unsigned int idx_a, unsigned int idx_b, unsigned int idx_c, unsigned int idx_d);

        std::map<unsigned int, std::vector<unsigned int> > m_bond_neigh_bonds; //!< Dictionary of bond index and four neighbor edges' indices
        std::map<unsigned int, unsigned int> m_bond_dihedral; //!< Dictionary of bond index and its dihedral index, 32767 if bond is on boundary
        std::map<unsigned int, std::vector<unsigned int>> m_bond_faces;
        std::vector<unsigned int> m_accepted_bonds;

        std::map<unsigned int, unsigned int> m_vertex_bond_num; //!< Dictionary of vertex index and number of connected bonds
        std::shared_ptr<BondData> m_bonds;
        std::shared_ptr<AngleData> m_angles;
        std::shared_ptr<DihedralData> m_dihedrals;
        // const std::shared_ptr<ComputeThermo> m_thermo;
        // const std::shared_ptr<Logger> m_logger; 
        // std::vector<std::string> m_loggable_quantities;
    };

namespace detail
    {
        //! Export the BondFlipUpdater class to python
        void export_BondFlipUpdater(pybind11::module& m);

    } // end namespace detail

// Third, this class offers a GPU accelerated method in order to demonstrate how to include CUDA
// code in pluins we need to declare a separate class for that (but only if ENABLE_HIP is set)

#ifdef ENABLE_HIP

//! A GPU accelerated nonsense particle updater written to demonstrate how to write a plugin w/ CUDA
//! code
/*! This updater simply sets all of the particle's velocities to 0 (on the GPU) when update() is
 * called.
 */
class BondFlipUpdaterGPU : public BondFlipUpdater
    {
    public:
    //! Constructor
    BondFlipUpdaterGPU(std::shared_ptr<SystemDefinition> sysdef, 
                    std::shared_ptr<Trigger> trigger,
                    // std::shared_ptr<ComputeThermo> thermo,
                    // std::shared_ptr<Logger> logger,
                    std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
                    std::map<unsigned int, unsigned int> bond_dihedral,
                    std::map<unsigned int, std::vector<unsigned int>> bond_faces);

    //! Take one timestep forward
    virtual void update(uint64_t timestep);
    };

namespace detail
    {
//! Export the BondFlipUpdaterGPU class to python
void export_BondFlipUpdaterGPU(pybind11::module& m);

    } // end namespace detail

#endif // ENABLE_HIP

    } // end namespace hoomd

#endif // _EXAMPLE_UPDATER_H_
