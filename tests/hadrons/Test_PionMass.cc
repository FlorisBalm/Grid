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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1,65536);
    std::stringstream seed_ss;
    seed_ss << dis(gen) << " " << dis(gen) << " " << dis(gen) << " " << dis(gen);
    std::string seed = seed_ss.str();

    // global parameters

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3425;
    globalPar.trajCounter.end   = 3430;
    globalPar.trajCounter.step  = 5;
    globalPar.seed              = seed;

    application.setPar(globalPar);
    // gauge field
//    MGauge::Load::Par loadPar;
 //   loadPar.file = "/home/floris/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Unit>("gauge");

    MSource::Z2::Par z2par;
    z2par.tA=0;
    z2par.tB=0;
    application.createModule<MSource::Z2>("Z2", z2par);

    MSource::Point::Par pointPar;
    pointPar.position="0 0 0 0";
    application.createModule<MSource::Point>("Pt",pointPar);

    MSource::U1::Par u1par;
    u1par.tA = 0;
    u1par.tB = 0;
    u1par.noise="548 6846 246 32";
    u1par.mom = "0 0 0 0";
    application.createModule<MSource::U1>("U1",u1par);

    time_t t = time(0);
    struct tm* now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%Y%m%d",now);
    std::string current_date(buffer);
    strftime(buffer, 80, "%H%M%S",now);
    std::string current_time(buffer); 

    for(double mass = -0.7; mass < .7; mass += 0.1){

        char buf[50];
        sprintf(buf, "%.4f", mass);
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
        quarkPar1.source = "U1";
        application.createModule<Quark>("QPt_"+mass_str, quarkPar1);


        Gamma g5(Gamma::Algebra::Gamma5);

        MContraction::Meson::Par mesPar;
        mesPar.q1 = "QPt_"+mass_str;
        mesPar.q2 = "QPt_"+mass_str;
        mesPar.gammaSource = Gamma::Algebra::Gamma5;
        mesPar.gammaSink = Gamma::Algebra::Gamma5;
        mesPar.output= "pionmass/"+current_date+"/PionMass_"+current_time+"_m"+mass_str;

        application.createModule<MContraction::Meson>("PionMass_"+mass_str,
                mesPar);
    }

    // execution
    application.saveParameterFile("pionmass/"+current_date+"/xml/PiMass"+current_time+".xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
