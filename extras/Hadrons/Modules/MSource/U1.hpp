/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSource/U1.hpp

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

#ifndef Hadrons_U1_hpp_
#define Hadrons_U1_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 U1 stochastic source
 -----------------------------
 At each point multiply with a stochastic U1 noise factor an inject momentum.
 Noise is uniform on unit circle
  
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - q: momentum (std::vector<RealD>)
 
 */
 
/******************************************************************************
 *                          U1 stochastic source                              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class U1Par: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(U1Par,
                                    unsigned int, tA,
                                    unsigned int, tB,
                                    std::string, mom);
};

template <typename FImpl>
class TU1: public Module<U1Par>
{
public:
    TYPE_ALIASES(FImpl,);
public:
    // constructor
    TU1(const std::string name);
    // destructor
    virtual ~TU1(void) = default;
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(U1, TU1<FIMPL>, MSource);

/******************************************************************************
 *                       TU1 template implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TU1<FImpl>::TU1(const std::string name)
: Module<U1Par>(name)
{
}


// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TU1<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TU1<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TU1<FImpl>::setup(void)
{
    env().template registerLattice<PropagatorField>(getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TU1<FImpl>::execute(void)
{
    Lattice<iScalar<vInteger>> t(env().getGrid());
    LatticeComplex eta (env().getGrid());

    LatticeComplex coord(env().getGrid());
    LatticeComplex qy(env().getGrid());
    qy=zero;
    std::vector<Real> momentum = strToVec<Real>(par().mom);
    for(unsigned int mu = 0; mu < Nd - 1; ++mu){
       LatticeCoordinate(coord,mu);
       qy += momentum[mu]*coord;
    }

    if (par().tA == par().tB)
    {
        LOG(Message) << "Generating U_1 wall source at t= " << par().tA
        << std::endl;
    }
    else
    {
        LOG(Message) << "Generating U_1 band for " << par().tA << " <= t <= "
        << par().tB << std::endl;
    }


    LatticeCoordinate(t,Tp);
    PropagatorField &src = *env().template createLattice<PropagatorField>(getName());
    random(*env().get4dRng(), eta);
    RealD TwoPi = (2.*M_PI);
    eta = eta*TwoPi;
    eta = where((t >= par().tA) and (t <= par().tB), eta, 0.*eta);
    eta = exp(timesI(eta));
    qy = exp(timesI(qy));
    src = 1.;
    src = src*qy*eta;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_U1_hpp_
