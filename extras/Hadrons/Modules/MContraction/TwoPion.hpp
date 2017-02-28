/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MContraction/TwoPion.hpp

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

#ifndef Hadrons_TwoPion_hpp_
#define Hadrons_TwoPion_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TTwoPion                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class TwoPionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TwoPionPar,
                                    std::string,    q1,
                                    std::string,    q2,
                                    std::string,    q3,
                                    std::string,    q4,
                                    std::string,    mom,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
class TTwoPion: public Module<TwoPionPar>
{
public:
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    TYPE_ALIASES(FImpl3, 3);
    TYPE_ALIASES(FImpl4, 4);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result, std::vector<Complex>, corr);
    };
public:
    // constructor
    TTwoPion(const std::string name);
    // destructor
    virtual ~TTwoPion(void) = default;
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_NS(TwoPion, ARG(TTwoPion<FIMPL,FIMPL,FIMPL,FIMPL>), MContraction);

/******************************************************************************
 *                           TTwoPion implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4>::TTwoPion(const std::string name)
: Module<TwoPionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().q3, par().q4};
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
std::vector<std::string> TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4>
void TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "' and '" 
                 << par().q3 << "' and '" << par().q4
                 << std::endl;
    XmlWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    PropagatorField2      &q2 = *env().template getObject<PropagatorField2>(par().q2);
    PropagatorField3      &q3 = *env().template getObject<PropagatorField3>(par().q3);
    PropagatorField4      &q4 = *env().template getObject<PropagatorField4>(par().q4);
    LatticeComplex        c(env().getGrid()),
                          d(env().getGrid());
    std::vector<TComplex> buf1,buf2;
    Result                result;

    LOG(Debug) << 126 << std::endl;
    LatticeComplex coord(env().getGrid());
    LatticeComplex py(env().getGrid());
    LatticeComplex pyneg(env().getGrid());
    
    py=zero;
    pyneg=zero;
    std::vector<Real> momentum = strToVec<Real>(par().mom);
    for(unsigned int mu = 0; mu < env().getNd() - 1; ++mu){
       LatticeCoordinate(coord,mu);
       py += momentum[mu]*coord;
    }
    
    pyneg = exp(-timesI(py));
    py = exp(timesI(py));
    LOG(Debug) << 137 << std::endl;
    
    c = trace(py*adj(q1)*q2);
    d = trace(pyneg*adj(q3)*q4);

    LOG(Debug) << 140 << std::endl; 
    sliceSum(c, buf1, Tp);
    sliceSum(d,buf2,Tp);
    result.corr.resize(buf1.size());
    
    for (unsigned int t = 0; t < buf1.size(); ++t)
    {
        result.corr[t] = TensorRemove(buf1[t])*TensorRemove(buf2[t]);
    }
    write(writer, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TwoPion_hpp_
