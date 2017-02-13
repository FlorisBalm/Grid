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
#include <Grid/Hadrons/Modules/MContraction/TwoPion.hpp>
#include <ctime>
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
    std::string flavour = "l";
    std::string down = "d";
    std::string up = "u";
    double mass = 1;
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1450;
    globalPar.trajCounter.end   = 1500;
    globalPar.trajCounter.step  = 10;
    globalPar.seed              = "0 2 3 4";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // sources
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);

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
    quarkPar.source = "pt";
    application.createModule<Quark>("Qpt_" + flavour, quarkPar);

    MContraction::TwoPion::Par mesPar;

    //Find current time for file to write to
    time_t t = time(0);
    struct tm* now = localtime(&t);
    char buffer[80];
    strftime(buffer, 80, "%Y%m%d-%H%M%S",now);
    std::string current_time(buffer);
    
    mesPar.output = "pions/pt_" + flavour + "_" + current_time; 
    mesPar.q1     = "Qpt_" + flavour;
    application.createModule<MContraction::TwoPion>("twopi_pt_"
                                                      + flavour,
                                                      mesPar);
    
    // execution
    application.saveParameterFile("TwoPion.xml");
    application.run();
     
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
