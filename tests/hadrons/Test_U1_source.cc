/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
 Copyright (C) 2015
 
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
 
 See the full license in the file "LICENSE" in the top level distribution
 directory.
 *******************************************************************************/

#include <Grid/Hadrons/Application.hpp>
#include <Grid/Hadrons/Modules/MSource/U1.hpp>
using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;
    
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;

    std::string flavour= "l";
    double mass = 0.1; 
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3425;
    globalPar.trajCounter.end   = 3445;
    globalPar.trajCounter.step  = 20;
    globalPar.seed              = "1 2 3 4";
    application.setPar(globalPar);
    // gauge field
    MGauge::Load::Par loadPar;
    loadPar.file = "/home/floris/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Load>("gauge", loadPar);
    
    std::vector<int> latt_size=GridDefaultLatt();

    RealD twoPiL = 2.*M_PI/double(latt_size[3]);
    MSource::U1::Par u1Par;
    u1Par.tA=5;
    u1Par.tB=5;
    
    u1Par.q={0,0,twoPiL};
    application.createModule<MSource::U1>("u1", u1Par);
    MAction::Wilson::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass  = mass;
    application.createModule<MAction::Wilson>("Wilson_" + flavour, actionPar);

    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action   = "Wilson_" + flavour;
    solverPar.residual = 1.0e-8;
    application.createModule<MSolver::RBPrecCG>("CG_" + flavour,
            solverPar);

    // propagators
    Quark::Par quarkPar;
    quarkPar.solver = "CG_" + flavour;
    quarkPar.source = "u1";
    application.createModule<Quark>("QU1_" + flavour, quarkPar);

    MContraction::Meson::Par mesPar;

    mesPar.output = "mesons/U1_" + flavour + flavour;
    mesPar.q1     = "QU1_" + flavour;
    mesPar.q2     = "QU1_" + flavour;
    application.createModule<MContraction::Meson>("meson_U1_"
            + flavour + flavour,
            mesPar);

    // execution
    application.saveParameterFile("U1.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
