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
#include <Grid/Hadrons/Modules/MSource/StochasticQuark.hpp>
#include <Grid/Hadrons/Modules/MContraction/RhoRho.hpp>
#include <Grid/Hadrons/Modules/MContraction/RhoPi.hpp>
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
    /*
    std::string in1,in2,in3,in4;
    std::cin >> in1,in2,in3,in4; 
    std::string seed = in1+"_"+in2+"_"+in3+"_"+in4; 
    */
    // global parameters

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3425;
    globalPar.trajCounter.end   = 3435;
    globalPar.trajCounter.step  = 5;
    globalPar.seed              = "1 2 432 125";

    application.setPar(globalPar);
    // gauge field
    MGauge::Load::Par loadPar;
    loadPar.file = "/home/floris/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Load>("gauge", loadPar);



    double mass = 0.1;
    char buf[50];
    sprintf(buf, "%.2f", mass);
    std::string mass_str(buf);

    //Wilson Action
    MAction::Wilson::Par actionPar;
    actionPar.gauge = "gauge";
    actionPar.mass  = mass;
    application.createModule<MAction::Wilson>("Wilson", actionPar);

    // solvers
    MSolver::RBPrecCG::Par solverPar;
    solverPar.action   = "Wilson";
    solverPar.residual = 1.0e-8;
    application.createModule<MSolver::RBPrecCG>("CG",solverPar);
    

    //Sources for quarks
    auto latt_size=GridDefaultLatt();
    RealD twoPiL = 2.*M_PI/double(latt_size[3]);
    std::vector<RealD> momentum = {0.,0.,twoPiL}; 
    std::string positive_momentum = "0 0 " + std::to_string(twoPiL);
    std::string negative_momentum = "0 0 " + std::to_string(-twoPiL);
    std::string zero_momentum = "0 0 0";
    
    std::string noise = "581 291 12 3";
    
    MSource::U1::Par sourcePar_pos;
    sourcePar_pos.tA = 0;
    sourcePar_pos.tB = 0;
    sourcePar_pos.mom= positive_momentum;
    sourcePar_pos.noise = noise;
    application.createModule<MSource::U1>("u1_p", sourcePar_pos);

    MSource::U1::Par sourcePar_neg;
    sourcePar_neg.tA = 0;
    sourcePar_neg.tB = 0;
    sourcePar_neg.mom= negative_momentum;
    sourcePar_neg.noise = noise;
    application.createModule<MSource::U1>("u1_n", sourcePar_neg);

    MSource::U1::Par sourcePar_zero;
    sourcePar_zero.tA = 0;
    sourcePar_zero.tB = 0;
    sourcePar_zero.mom= zero_momentum;
    sourcePar_zero.noise = noise;
    application.createModule<MSource::U1>("u1_0", sourcePar_zero);


    Quark::Par quarkPos;
    quarkPos.source = "u1_p";
    quarkPos.solver = "CG";
    application.createModule<Quark>("Qu1_p", quarkPos);

    Quark::Par quarkNeg;
    quarkNeg.source = "u1_n";
    quarkNeg.solver = "CG";
    application.createModule<Quark>("Qu1_n", quarkNeg);

    Quark::Par quarkZero;
    quarkZero.source = "u1_0";
    quarkZero.solver="CG";
    application.createModule<Quark>("Qu1_0", quarkZero);

    MSource::StochasticQuark::Par wtype1par;
    wtype1par.q = "Qu1_p";
    wtype1par.tA = 0;
    wtype1par.tB = 0;
    wtype1par.mom = negative_momentum;
    application.createModule<MSource::StochasticQuark>("Stoc1", wtype1par);

    MSource::StochasticQuark::Par wtype2par;
    wtype2par.q = "Qu1_n";
    wtype2par.tA = 0;
    wtype2par.tB = 0;
    wtype2par.mom = negative_momentum;
    application.createModule<MSource::StochasticQuark>("Stoc2", wtype1par);

    // propagators
     
    Quark::Par quarkPar_S1;
    quarkPar_S1.solver = "CG";
    quarkPar_S1.source = "Stoc1";
    application.createModule<Quark>("QS1", quarkPar_S1);

    Quark::Par quarkPar_S2;
    quarkPar_S2.solver = "CG";
    quarkPar_S2.source = "Stoc2";
    application.createModule<Quark>("QS2", quarkPar_S2);

    time_t t = time(0);
    struct tm* now = localtime(&t);

    char buffer[80];
    strftime(buffer, 80, "%Y%m%d",now);
    std::string current_date(buffer);
    strftime(buffer, 80, "%H%M%S",now);
    std::string current_time(buffer); 

    MContraction::RhoPi::Par rhoPiPar;
    rhoPiPar.output= "rhopi/"+current_date+"/RP_"+current_time;
    rhoPiPar.q1 = "Qu1_0";
    rhoPiPar.q2 = "QS1";
    rhoPiPar.q3 = "QS2";
    application.createModule<MContraction::RhoPi>("rhopi",rhoPiPar);

    // execution
    application.saveParameterFile("RhoPi.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
