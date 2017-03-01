/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSource/StochasticQuark.hpp

Copyright (C) 2015
Copyright (C) 2016

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_StochasticQuark_hpp_
#define Hadrons_StochasticQuark_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Sequential source
 -----------------------------
 * src_x = q_x * theta(x_3 - tA) * theta(tB - x_3) * gamma * exp(i x.mom)
 
 * options:
 - q: input propagator (string)
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - gamma: gamma product to insert (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         StochasticQuark                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class StochasticQuarkPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StochasticQuarkPar,
                                    std::string,    q,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    std::string,    mom);
};

template <typename FImpl>
class TStochasticQuark: public Module<StochasticQuarkPar>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TStochasticQuark(const std::string name);
    // destructor
    virtual ~TStochasticQuark(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(StochasticQuark, TStochasticQuark<FIMPL>, MSource);

/******************************************************************************
 *                         TStochasticQuark implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStochasticQuark<FImpl>::TStochasticQuark(const std::string name)
: Module<StochasticQuarkPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStochasticQuark<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStochasticQuark<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStochasticQuark<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStochasticQuark<FImpl>::execute(void)
{
    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source at t= " << par().tA << std::endl;
    }
    else
    {
        LOG(Message) << "Generating gamma_" << par().gamma
                     << " sequential source for "
                     << par().tA << " <= t <= " << par().tB << std::endl;
    }
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    PropagatorField &q   = *env().template getObject<PropagatorField>(par().q);

    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex             ph(env().getGrid()),
                               coor(env().getGrid());
    Gamma                      g5(Algebra::Gamma5);
    std::vector<Real>          p;
    Complex                    i(0.0,1.0);
    
    p  = strToVec<Real>(par().mom);
    ph = zero;
    for(unsigned int mu = 0; mu < env().getNd()-1; mu++)
    {
        LatticeCoordinate(coor, mu);
        ph += p[mu]*coor;
    }
    ph = exp(i*ph);
    LatticeCoordinate(t, Tp);
    src = where((t >= par().tA) and (t <= par().tB), ph*(g5*q), 0.*q);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_StochasticQuark_hpp_
