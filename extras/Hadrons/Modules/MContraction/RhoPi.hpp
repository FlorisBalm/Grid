/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/RhoPi.hpp

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

#ifndef Hadrons_RhoPi_hpp_
#define Hadrons_RhoPi_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TRhoPi                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class RhoPiPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RhoPiPar,
                                    std::string,    q1,
                                    std::string,    q2,
                                    std::string,    q3,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3>
class TRhoPi: public Module<RhoPiPar>
{
public:
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    TYPE_ALIASES(FImpl3, 3);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result, std::vector<std::vector<Complex>>, corr);
    };
public:
    // constructor
    TRhoPi(const std::string name);
    // destructor
    virtual ~TRhoPi(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(RhoPi, ARG(TRhoPi<FIMPL, FIMPL,FIMPL>), MContraction);

/******************************************************************************
 *                           TRhoPi implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
TRhoPi<FImpl1, FImpl2,FImpl3>::TRhoPi(const std::string name)
: Module<RhoPiPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TRhoPi<FImpl1, FImpl2,FImpl3>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3>
std::vector<std::string> TRhoPi<FImpl1, FImpl2,FImpl3>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3>
void TRhoPi<FImpl1, FImpl2,FImpl3>::execute(void)
{
    LOG(Message) << "Computing pipi-rho contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "' and '" << par().q3
                 << std::endl;
    
    XmlWriter             writer(par().output);

    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    PropagatorField3      &q3 = *env().template getObject<PropagatorField3>(par().q3);

    LatticeComplex        c1(env().getGrid()),
                          c2(env().getGrid());

    Gamma                 g3(Gamma::Algebra::GammaZ);
    Gamma                 g5(Gamma::Algebra::Gamma5);

    std::vector<TComplex> buf1, buf2;
    Result                result;
    
    c1 = trace(q1*g5*g3*adj(q2));
    c2 = trace(q1*g5*g3*adj(q3));
    sliceSum(c1, buf1, Tp);
    sliceSum(c2, buf2, Tp);


    result.corr.push_back(std::vector<Complex>(buf1.size()));
    result.corr.push_back(std::vector<Complex>(buf1.size()));
    result.corr.push_back(std::vector<Complex>(buf1.size()));
    
    for (unsigned int t = 0; t < buf1.size(); ++t)
    { 
        result.corr[0][t] = TensorRemove(buf1[t]);
        result.corr[1][t] = TensorRemove(buf2[t]);
        result.corr[2][t] = TensorRemove(buf1[t]) + TensorRemove(buf2[t]);
    }
    write(writer, "rhopi", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_RhoPi_hpp_
