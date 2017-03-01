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
                                    std::string,    q_pos,
                                    std::string,    q_neg,
                                    std::string,    q0_1,
                                    std::string,    q0_2,
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
    std::vector<std::string> input = {par().q_pos, par().q_neg, par().q0_1, par().q0_2};
    
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
                 << " quarks '" << par().q_pos << "' and '" << par().q_neg << "' and '" 
                 << par().q0_1 << "' and '" << par().q0_2
                 << std::endl;
    XmlWriter             writer(par().output);
    PropagatorField1      &q_pos = *env().template getObject<PropagatorField1>(par().q_pos);
    PropagatorField2      &q_neg = *env().template getObject<PropagatorField2>(par().q_neg);
    PropagatorField3      &q0_1 = *env().template getObject<PropagatorField3>(par().q0_1);
    PropagatorField4      &q0_2 = *env().template getObject<PropagatorField4>(par().q0_2);
    LatticeComplex        traceCrossed1(env().getGrid()),
                          traceCrossed2(env().getGrid()),
                          tracePaired1(env().getGrid()),
                          tracePaired2(env().getGrid());
    std::vector<std::vector<TComplex>> buf;
    Result                result;

    LatticeComplex coord(env().getGrid());
    LatticeComplex py(env().getGrid());
    LatticeComplex pyneg(env().getGrid());
    
    py=zero;
    pyneg=zero;
    std::vector<Real> momentum = strToVec<Real>(par().mom);
    for(unsigned int mu = 0; mu < env().getNd() - 1; ++mu){
       LatticeCoordinate(coord,mu);
       LOG(Debug) << coord << std::endl;
       py += momentum[mu]*coord;
    }
    
    pyneg = exp(-timesI(py));
    py = exp(timesI(py));
    LOG(Debug) << py << "\n" << pyneg << std::endl;
    //Paired diagrams diagram 
    for(int i = 0; i < 4; ++i){
        buf.push_back(std::vector<TComplex>());
    }

    traceCrossed1 = trace(py*adj(q0_1)*q_neg);
    traceCrossed2 = trace(pyneg*adj(q0_2)*q_pos);
    
    sliceSum(traceCrossed1,buf[0],Tp);
    sliceSum(traceCrossed2,buf[1],Tp);

    tracePaired1 = trace(py*adj(q0_1)*q_pos);
    tracePaired2 = trace(pyneg*adj(q0_2)*q_neg);

    sliceSum(tracePaired1,buf[2],Tp);
    sliceSum(tracePaired2,buf[3],Tp);


    //Other traces
    result.corr.resize(buf[0].size());
    
    LOG(Message) << " Got here to where I put stuff in" << std::endl;
    for (unsigned int t = 0; t < buf[0].size(); ++t)
    {
        auto tot = TensorRemove(buf[0][t]); 
        LOG(Message) << " Still no crash" << std::endl;
        for(int i = 1; i < 4 ; ++i){
            LOG(Message) << " Still not (2)" << std::endl;
            if(buf[i].size() != 0){
                LOG(Message) << " Not (3)" << std::endl;
                tot *= TensorRemove(buf[i][t]);
            }
        }
        result.corr[t] = tot;
    }
    write(writer, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TwoPion_hpp_
