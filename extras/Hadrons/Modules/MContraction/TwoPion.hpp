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
                                std::string, q1,
                                std::string, output);
};

template <typename FImpl1>
class TTwoPion: public Module<TwoPionPar>
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
                TTwoPion(const std::string name);
                // destructor
                virtual ~TTwoPion(void) = default;
                // dependencies/products
                virtual std::vector<std::string> getInput(void);
                virtual std::vector<std::string> getOutput(void);
                // execution
                virtual void execute(void);
};

MODULE_REGISTER_NS(TwoPion, ARG(TTwoPion<FIMPL>), MContraction);

/******************************************************************************
 *                           TTwoPion implementation                            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
        template <typename FImpl1>
TTwoPion<FImpl1>::TTwoPion(const std::string name)
        : Module<TwoPionPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
        template <typename FImpl1>
std::vector<std::string> TTwoPion<FImpl1>::getInput(void)
{
        std::vector<std::string> input = {par().q1};

        return input;
}

        template <typename FImpl1>
std::vector<std::string> TTwoPion<FImpl1>::getOutput(void)
{
        std::vector<std::string> output = {getName()};


        return output;
}
// execution ///////////////////////////////////////////////////////////////////
        template <typename FImpl1>
void TTwoPion<FImpl1>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << std::endl;
    
    XmlWriter             writer(par().output);
    PropagatorField1      &q1 = *env().template getObject<PropagatorField1>(par().q1);
    LatticeComplex        c(env().getGrid());
    Gamma                 g5(Gamma::Algebra::Gamma5);
    std::vector<TComplex> buf;
    Result                result;

    c = (trace(q1*g5*g5*adj(q1)*g5*g5)*trace(q1*g5*g5*adj(q1)*g5*g5) + trace(q1*g5*g5*adj(q1)*g5*g5*q1*g5*g5*adj(q1)*g5*g5))+(trace(q1*g5*g5*adj(q1)*g5*g5)*trace(q1*g5*g5*adj(q1)*g5*g5) + trace(q1*g5*g5*adj(q1)*g5*g5*q1*g5*g5*adj(q1)*g5*g5));
    sliceSum(c, buf, Tp);
    result.corr.resize(buf.size());
    for (unsigned int t = 0; t < buf.size(); ++t)
    {
            result.corr[t] = TensorRemove(buf[t]);
    }
    write(writer, "twopion", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_TwoPion_hpp_
