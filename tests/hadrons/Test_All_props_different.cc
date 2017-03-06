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
#include <Grid/Hadrons/Modules/MContraction/TwoPion.hpp>
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

    // global parameters

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3425;
    globalPar.trajCounter.end   = 3430;
    globalPar.trajCounter.step  = 5;
    globalPar.seed              = "1 4 3 2";
    
    application.setPar(globalPar);
    // gauge field
    MGauge::Load::Par loadPar;
    loadPar.file = "/home/floris/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Load>("gauge", loadPar);
    

    double mass = 0.1;
    
    auto latt_size=GridDefaultLatt();
    RealD twoPiL = 2.*M_PI/double(latt_size[3]);
    std::string momentum = "0 0 " + std::to_string(twoPiL);
    std::string negative_momentum = "0 0 " + std::to_string(-twoPiL);

    for(int i = 1; i<3; ++i){
        std::string n = std::to_string(i);
        
        MSource::U1::Par u1Par_posMomentum;
        u1Par_posMomentum.tA=0;
        u1Par_posMomentum.tB=0;
        u1Par_posMomentum.mom = momentum;
        application.createModule<MSource::U1>("u1_p_"+n, u1Par_posMomentum);


        MSource::U1::Par u1Par_negMomentum;
        u1Par_negMomentum.tA=0;
        u1Par_negMomentum.tB=0;
        u1Par_negMomentum.mom = negative_momentum;

        application.createModule<MSource::U1>("u1_n_"+n, u1Par_negMomentum);


        MSource::U1::Par u1Par_zeroMomentum;
        u1Par_zeroMomentum.tA=0;
        u1Par_zeroMomentum.tB=0;
        u1Par_zeroMomentum.mom = "0 0 0";

        application.createModule<MSource::U1>("u1_0_"+n, u1Par_zeroMomentum);

        MSource::U1::Par u1Par_zeroMomentum2;
        u1Par_zeroMomentum2.tA=0;
        u1Par_zeroMomentum2.tB=0;
        u1Par_zeroMomentum2.mom = "0 0 0";

        application.createModule<MSource::U1>("u1_0-2_"+n, u1Par_zeroMomentum2);
    } 

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

    std::vector<std::string> propagatorNames = {"u1_p", "u1_n","u1_0","u1_0-2"};
    for(auto s : propagatorNames){
        
        for(int i = 1; i<3; ++i){
            std::string n = std::to_string(i); 
            Quark::Par quarkPar;
            quarkPar.solver = "CG";
            quarkPar.source = s + "_" + n;
            application.createModule<Quark>("Q"+s+"_"+n, quarkPar);
        }

    }

    // propagators
    /*
    Quark::Par quarkPar1;
    quarkPar1.solver = "CG";
    quarkPar1.source = "u1_p";
    application.createModule<Quark>("QU1_p", quarkPar1);

    Quark::Par quarkPar2;
    quarkPar2.solver = "CG";
    quarkPar2.source = "u1_0";
    application.createModule<Quark>("QU1_0", quarkPar2);


    Quark::Par quarkPar3;
    quarkPar3.solver = "CG";
    quarkPar3.source = "u1_n";
    application.createModule<Quark>("QU1_n", quarkPar3);

    Quark::Par quarkPar4;
    quarkPar4.solver = "CG";
    quarkPar4.source = "u1_0";
    application.createModule<Quark>("QU1_0_2", quarkPar4);

    */
    
    time_t t = time(0);
    struct tm* now = localtime(&t);
    
    char buffer[80];
    strftime(buffer, 80, "%Y%m%d",now);
    std::string current_date(buffer);
    strftime(buffer, 80, "%H%M%S",now);
    std::string current_time(buffer); 

    MContraction::TwoPion::Par twoPionPar;
    twoPionPar.output = "twopion/"+current_date+"/U1_" + current_time;
    std::cout << twoPionPar.output << std::endl;
     
    twoPionPar.q0_1_1     = "Qu1_0_1";
    twoPionPar.q_pos_1     = "Qu1_p_1";
    twoPionPar.q0_2_1    = "Qu1_0-2_1";
    twoPionPar.q_neg_1     = "Qu1_n_1";
    twoPionPar.q0_1_2     = "Qu1_0_2";
    twoPionPar.q_pos_2     = "Qu1_p_2";
    twoPionPar.q0_2_2    = "Qu1_0-2_2";
    twoPionPar.q_neg_2     = "Qu1_n_2";

    twoPionPar.mom    = momentum;

    application.createModule<MContraction::TwoPion>("twoPion_U1",
            twoPionPar);

    // execution
    application.saveParameterFile("U1.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
