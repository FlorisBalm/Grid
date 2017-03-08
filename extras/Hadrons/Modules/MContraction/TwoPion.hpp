/*************************************************************************************
       py += momentum[mu]*coord;

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
                                    std::string,    q_pos_1,
                                    std::string,    q_neg_1,
                                    std::string,    q0_1_1,
                                    std::string,    q0_2_1,
                                    std::string,    q_pos_2,
                                    std::string,    q_neg_2,
                                    std::string,    q0_1_2,
                                    std::string,    q0_2_2,
                                    std::string,    qs_pn,
                                    std::string,    qs_np,
                                    std::vector<std::string>, v_qs_p0,
                                    std::string,    mom,
                                    std::string,    output);
};

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4,typename FImpl5, typename FImpl6, typename FImpl7, typename FImpl8, typename FImpl9, typename FImpl10, typename FImpl11, typename FImpl12>
class TTwoPion: public Module<TwoPionPar>
{
public:
    TYPE_ALIASES(FImpl1, 1);
    TYPE_ALIASES(FImpl2, 2);
    TYPE_ALIASES(FImpl3, 3);
    TYPE_ALIASES(FImpl4, 4);
    TYPE_ALIASES(FImpl5, 5);
    TYPE_ALIASES(FImpl6, 6);
    TYPE_ALIASES(FImpl7, 7);
    TYPE_ALIASES(FImpl8, 8);
    TYPE_ALIASES(FImpl9, 9);
    TYPE_ALIASES(FImpl10, 10);
    TYPE_ALIASES(FImpl11, 11);
    TYPE_ALIASES(FImpl12, 12);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result, std::vector<std::vector<Complex>>, corr);
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

MODULE_REGISTER_NS(TwoPion, ARG(TTwoPion<FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL,FIMPL>), MContraction);

/******************************************************************************
 *                           TTwoPion implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4,typename FImpl5, typename FImpl6, typename FImpl7, typename FImpl8, typename FImpl9, typename FImpl10, typename FImpl11, typename FImpl12>
TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4,FImpl5,FImpl6,FImpl7,FImpl8,FImpl9,FImpl10,FImpl11,FImpl12>::TTwoPion(const std::string name)
: Module<TwoPionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4,typename FImpl5, typename FImpl6, typename FImpl7, typename FImpl8, typename FImpl9, typename FImpl10, typename FImpl11, typename FImpl12>
std::vector<std::string> TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4,FImpl5,FImpl6,FImpl7,FImpl8,FImpl9,FImpl10,FImpl11,FImpl12>::getInput(void)
{
    std::vector<std::string> input = {par().q_pos_1, par().q_neg_1, par().q0_1_1, par().q0_2_1,par().q_pos_2, par().q_neg_2, par().q0_1_2, par().q0_2_2}
    
    return input;
}

template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4,typename FImpl5, typename FImpl6, typename FImpl7, typename FImpl8, typename FImpl9, typename FImpl10, typename FImpl11, typename FImpl12>
std::vector<std::string> TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4,FImpl5,FImpl6,FImpl7,FImpl8,FImpl9,FImpl10,FImpl11,FImpl12>::getOutput(void)
{
    std::vector<std::string> output = {getName()};
    
    return output;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2, typename FImpl3, typename FImpl4,typename FImpl5, typename FImpl6, typename FImpl7, typename FImpl8, typename FImpl9, typename FImpl10, typename FImpl11, typename FImpl12>
void TTwoPion<FImpl1,FImpl2,FImpl3,FImpl4,FImpl5,FImpl6,FImpl7,FImpl8,FImpl9,FImpl10,FImpl11,FImpl12>::execute(void)
{
//    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
 //j                << " quarks '" << par().q_pos << "' and '" << par().q_neg << "' and '" 
   //              << par().q0_1 << "' and '" << par().q0_2
    //             << std::endl;
    XmlWriter             writer(par().output);
    PropagatorField1      &q_pos_1 = *env().template getObject<PropagatorField1>(par().q_pos_1);
    PropagatorField2      &q_neg_1 = *env().template getObject<PropagatorField2>(par().q_neg_1);
    PropagatorField3      &q0_1_1 = *env().template getObject<PropagatorField3>(par().q0_1_1);
    PropagatorField4      &q0_2_1 = *env().template getObject<PropagatorField4>(par().q0_2_1);
    PropagatorField5      &q_pos_2 = *env().template getObject<PropagatorField5>(par().q_pos_2);
    PropagatorField6      &q_neg_2 = *env().template getObject<PropagatorField6>(par().q_neg_2);
    PropagatorField7      &q0_1_2 = *env().template getObject<PropagatorField7>(par().q0_1_2);
    PropagatorField8      &q0_2_2 = *env().template getObject<PropagatorField8>(par().q0_2_2);
    PropagatorField9      &qs_np = *env().template getObject<PropagatorField9>(par().qs_np);
    PropagatorField10      &qs_pn = *env().template getObject<PropagatorField10>(par().qs_pn);



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
    for(int i = 0; i < 6; ++i){
        buf.push_back(std::vector<TComplex>());
    }

    traceCrossed1 = trace(py*adj(q0_1_1)*q_neg_1);
    traceCrossed2 = trace(pyneg*adj(q0_2_2)*q_pos_2);
    
    sliceSum(traceCrossed1,buf[0],Tp);
    sliceSum(traceCrossed2,buf[1],Tp);

    tracePaired1 = trace(py*adj(q0_1_2)*q_pos_2);
    tracePaired2 = trace(pyneg*adj(q0_2_2)*q_neg_2);

    sliceSum(tracePaired1,buf[2],Tp);
    sliceSum(tracePaired2,buf[3],Tp);

          
    Lattice<iScalar<vInteger>> timecoords(env().getGrid());
    LatticeCoordinate(timecoords,Tp);
    PropagatorField11 &qs_p0 = *env().template getObject<PropagatorField11>(par().v_qs_p0[0]);
    qs_p0 = where((timecoords == 0), qs_p0, 0*qs_p0);
    for (int t = 1; t < Td; ++t){
        std::string qname = par().v_qs_p0[t];

        PropagatorField12 &qs_p0_t =  *env().template getObject<PropagatorField12>(par().v_qs_p0[t]);
        qs_p0_t = where(t==timecoords, qs_p0_t, 0*qs_p0_t);
        qs_p0 += qs_p0_t;
    }
    traceSquare = trace(adj(qs_pn)*pyneg*qs_p0);
    traceHourglass = trace(adj(qs_np)pyneg*qs_p0);
    
    sliceSum(traceSquare,buf[4],Tp);
    sliceSum(traceHourglass,buf[5],Tp);

    //Other traces
    result.corr.push_back(std::vector<Complex>(buf[0].size()));
    result.corr.push_back(std::vector<Complex>(buf[0].size()));
    result.corr.push_back(std::vector<Complex>(buf[0].size()));
    result.corr.push_back(std::vector<Complex>(buf[0].size()));
    result.corr.push_back(std::vector<Complex>(buf[0].size()));
    for (unsigned int t = 0; t < buf[0].size(); ++t)
    {
        result.corr[0][t] = -TensorRemove(buf[2][t])*TensorRemove(buf[3][t]); 
        result.corr[1][t] =  TensorRemove(buf[0][t])*TensorRemove(buf[1][t]);
        result.corr[2][t] =  TensorRemove(buf[4][t]);
        result.corr[3][t] = -TensorRemove(buf[5][t]);
        result.corr[4][t] = 2*result.corr[3][t] + 2*result.corr[2][t] + result.corr[1][t] + result.corr[0][t];
    }
    write(writer, "twopion", result);
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TwoPion_hpp_
