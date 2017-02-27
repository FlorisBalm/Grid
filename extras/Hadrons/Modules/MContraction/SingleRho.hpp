/*************************************************************************************

  Grid physics library, www.github.com/paboyle/Grid 

  Source file: extras/Hadrons/Modules/MContraction/Meson.hpp

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

#ifndef Hadrons_SingleRho_hpp_
#define Hadrons_SingleRho_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp> 
BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                TSingleRho                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)
        class SingleRhoPar: Serializable
{
        public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(SingleRhoPar,
                                std::string, q1,
                                std::string, q2,
                                std::string, q3,
                                std::string, q4,
                                std::string, output);
};

template <typename FImpl1>
class TSingleRho: public Module<SingleRhoPar>
{
        public:
                TYPE_ALIASES(FImpl1, 1);

                class Result: Serializable
        {
                public:
                        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr);
        };
        public:
                // constructor
                TSingleRho(const std::string name);
                // destructor
                virtual ~TSingleRho(void) = default;
                // dependencies/products
                virtual std::vector<std::string> getInput(void);
                virtual std::vector<std::string> getOutput(void);
                // execution
                virtual void execute(void);
};

MODULE_REGISTER_NS(SingleRho, ARG(TSingleRho<FIMPL>), MContraction);

/******************************************************************************
 *                           TSingleRho implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1>
TSingleRho< FImpl1>::TSingleRho(const std::string name)
        : Module<SingleRhoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1>
std::vector<std::string> TSingleRho< FImpl1>::getInput(void)
{
        std::vector<std::string> input = {par().q1};

        return input;
}

template <typename FImpl1>
std::vector<std::string> TSingleRho< FImpl1>::getOutput(void)
{
        std::vector<std::string> output = {getName()};


        return output;
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1>
void TSingleRho< FImpl1>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '"
                 << par().q1 << ", "  
                 << std::endl;
    
    XmlWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    LatticeComplex        c(env().getGrid());
    SpinMatrix            g5, g3;
    std::vector<TComplex> buf;
    Result                result;
    g3 = makeGammaProd(1<<3); 
    g5 = makeGammaProd(Ns*Ns - 1);
    c=trace(q1*g3*q1*g3);
    //c = (trace(q1*g5*g5*adj(q1)*g5*g5)*trace(q1*g5*g5*adj(q1)*g5*g5) + trace(q1*g5*g5*adj(q1)*g5*g5*q1*g5*g5*adj(q1)*g5*g5))+(trace(q1*g5*g5*adj(q1)*g5*g5)*trace(q1*g5*g5*adj(q1)*g5*g5) + trace(q1*g5*g5*adj(q1)*g5*g5*q1*g5*g5*adj(q1)*g5*g5));
    sliceSum(c, buf, Tp);
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
            result.corr[t] = TensorRemove(buf[t]);
    }


    write(writer, "singlerho", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_SingleRho_hpp_
