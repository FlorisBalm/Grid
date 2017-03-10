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
#include <Grid/Hadrons/Modules/MSource/Z2.hpp>
#include <Grid/Hadrons/Modules/MContraction/RhoRho.hpp>
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
    std::stringstream seed_ss;
    seed_ss << rand() << " " << rand() << " " << rand() << " " << rand();
    std::string seed = seed_ss.str();

    // global parameters

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3425;
    globalPar.trajCounter.end   = 3435;
    globalPar.trajCounter.step  = 5;
    globalPar.seed              = seed;

    application.setPar(globalPar);
    // gauge field
    MGauge::Load::Par loadPar;
    loadPar.file = "/home/floris/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Load>("gauge", loadPar);
    MSource::Z2::Par z2par;
    z2par.tA=0;
    z2par.tB=0;
    application.createModule<MSource::Z2>("Z2", z2par);

    double mass = 0.1;
    char buf[50];
    sprintf(buf, "%.2f", mass);
    std::string mass_str(buf);

    //Wilson Action
    MAction::Wilson::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass  = mass;
    application.createModule<MAction::Wilson>("Wilson_" + mass_str, actionPar);

    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action   = "Wilson_"+mass_str;
    solverPar.residual = 1.0e-8;
    application.createModule<MSolver::RBPrecCG>("CG_"+mass_str,solverPar);

    // propagators
    Quark::Par quarkPar1;
    quarkPar1.solver = "CG_"+mass_str;
    quarkPar1.source = "Z2";
    application.createModule<Quark>("QZ2_"+mass_str, quarkPar1);


    time_t t = time(0);
    struct tm* now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%Y%m%d",now);
    std::string current_date(buffer);
    strftime(buffer, 80, "%H%M%S",now);
    std::string current_time(buffer); 

    MContraction::RhoRho::Par rhoPar;
    rhoPar.q1 = "QZ2_"+mass_str;
    rhoPar.q2 = "QZ2_"+mass_str;
    rhoPar.output= "rhorho/"+current_date+"/Z2_"+current_time+"_m"+mass_str+"_"+seed;

    application.createModule<MContraction::RhoRho>("RhoRho_Z2_"+mass_str,
            rhoPar);


    // execution
    application.saveParameterFile("rhorho/"+current_date+"/xml/RhoRho"+current_time+".xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
