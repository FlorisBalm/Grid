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
#include <Grid/Hadrons/Modules/MContraction/TwoPion.hpp>
#include <random>
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

    //Set up seeds. These are recorded in xml files 
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1,65536);
    
    std::stringstream seed1_ss, noise1_ss, noise2_ss;

    seed1_ss << dis(gen) << " " << dis(gen) << " " << dis(gen) << " " << dis(gen);
    noise1_ss << dis(gen) << " " << dis(gen) << " " << dis(gen) << " " << dis(gen);
    noise2_ss << dis(gen) << " " << dis(gen) << " " << dis(gen) << " " << dis(gen);

    std::string seed = seed1_ss.str(),
                noise1 = noise1_ss.str(),
                noise2 = noise2_ss.str();
    // run setup ///////////////////////////////////////////////////////////////
    Application              application;

    // global parameters

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 3445;
    globalPar.trajCounter.end   = 3450;
    globalPar.trajCounter.step  = 5;
    globalPar.seed              = seed;
    application.setPar(globalPar);
    // gauge field
    MGauge::Load::Par loadPar;
    loadPar.file = "/home/s1205916/mphys/configurations/ckpoint_lat";
    application.createModule<MGauge::Load>("gauge", loadPar);
    
    double mass = -0.769;
    
    auto latt_size=GridDefaultLatt();
    RealD twoPiL = 2.*M_PI/double(latt_size[Zp]);
    std::vector<RealD> momentum = {0.,0.,twoPiL}; 
    std::string positive_momentum = "0 0 " + std::to_string(twoPiL);
    std::string negative_momentum = "0 0 " + std::to_string(-twoPiL);
    std::string zero_momentum = "0 0 0";
    /*** ALTERNATIVE, BUT I THINK THERE"S NO ENV HERE YET...
    LatticeComplex noise1(env().getGrid()),
                   noise2(env().getGrid());
    random(*env().get4dRng(),noise1);
    random(*env().get4dRng(),noise2);

    noise1 = noise1* twoPi;
    noise2 = noise2* twoPi;
    noise1 = exp(timesI(noise1));
    noise2 = exp(timesI(noise2));

    ***/

    MSource::U1::Par u1Par_posMomentum_1;
    u1Par_posMomentum_1.tA=0;
    u1Par_posMomentum_1.tB=0;
    u1Par_posMomentum_1.mom = positive_momentum;
    u1Par_posMomentum_1.noise = noise1;
    application.createModule<MSource::U1>("u1_p_1", u1Par_posMomentum_1);


    MSource::U1::Par u1Par_negMomentum_1;
    u1Par_negMomentum_1.tA=0;
    u1Par_negMomentum_1.tB=0;
    u1Par_negMomentum_1.mom = negative_momentum;
    u1Par_negMomentum_1.noise = noise1; 
    application.createModule<MSource::U1>("u1_n_1", u1Par_negMomentum_1);


    MSource::U1::Par u1Par_zeroMomentum;
    u1Par_zeroMomentum.tA=0;
    u1Par_zeroMomentum.tB=0;
    u1Par_zeroMomentum.mom = zero_momentum;
    u1Par_zeroMomentum.noise=noise1;
    application.createModule<MSource::U1>("u1_0-1", u1Par_zeroMomentum);

    MSource::U1::Par u1Par_zeroMomentum2;
    u1Par_zeroMomentum2.tA=0;
    u1Par_zeroMomentum2.tB=0;
    u1Par_zeroMomentum2.mom = zero_momentum;
    u1Par_zeroMomentum2.noise=noise2;
    application.createModule<MSource::U1>("u1_0-2", u1Par_zeroMomentum2);

    MSource::U1::Par u1Par_posMomentum_2;
    u1Par_posMomentum_2.tA=0;
    u1Par_posMomentum_2.tB=0;
    u1Par_posMomentum_2.mom = positive_momentum;
    u1Par_posMomentum_2.noise = noise2;
    application.createModule<MSource::U1>("u1_p_2", u1Par_posMomentum_2);


    MSource::U1::Par u1Par_negMomentum_2;
    u1Par_negMomentum_2.tA=0;
    u1Par_negMomentum_2.tB=0;
    u1Par_negMomentum_2.mom = negative_momentum;
    u1Par_negMomentum_2.noise = noise2;
    application.createModule<MSource::U1>("u1_n_2", u1Par_negMomentum_2);



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

    std::vector<std::string> propagatorNames = {"u1_p_1","u1_p_2", "u1_n_1",  "u1_n_2", "u1_0-1","u1_0-2"};
    for(auto s : propagatorNames){
        Quark::Par quarkPar;
        quarkPar.solver = "CG";
        quarkPar.source = s;
        application.createModule<Quark>("Q"+s, quarkPar);

    }


    MSource::StochasticQuark::Par stoch_pos_neg;
    stoch_pos_neg.q = "Qu1_n_1"; //this needs to be 2 to get constistent noise, it's this now due to 
    stoch_pos_neg.tA = 0;
    stoch_pos_neg.tB = 0;
    stoch_pos_neg.mom = positive_momentum;
    application.createModule<MSource::StochasticQuark>("S_PN",stoch_pos_neg);
     
    MSource::StochasticQuark::Par stoch_neg_pos;
    stoch_neg_pos.q = "Qu1_p_1";
    stoch_neg_pos.tA = 0;
    stoch_neg_pos.tB = 0;
    stoch_neg_pos.mom = negative_momentum;
    application.createModule<MSource::StochasticQuark>("S_NP",stoch_neg_pos);
    
    Quark::Par sQuarkPar;
    sQuarkPar.solver="CG";
    sQuarkPar.source="S_PN";
    application.createModule<Quark>("QS_PN",sQuarkPar);
    
    Quark::Par sQuarkPar2;
    sQuarkPar2.solver="CG";
    sQuarkPar2.source="S_NP";
    application.createModule<Quark>("QS_NP",sQuarkPar2);

    std::vector<std::string> wQuarkNames;
    std::vector<std::string> wSourceNames;
    unsigned int Time = latt_size[Tp];
    for(unsigned int t = 0; t < Time; ++t){
        MSource::StochasticQuark::Par stoch_p_0;
        stoch_p_0.q = "Qu1_0-1";
        stoch_p_0.tA = t;
        stoch_p_0.tB = t;
        stoch_p_0.mom = zero_momentum;
        std::string stochname = "S_P0_"+std::to_string(t);
        application.createModule<MSource::StochasticQuark>(stochname,stoch_p_0);
        wSourceNames.push_back(stochname); 
    }
    for(unsigned int t = 0; t<wSourceNames.size(); ++t){
        Quark::Par quarkPar;
        quarkPar.solver="CG";

        quarkPar.source=wSourceNames[t];
        wQuarkNames.push_back("Q"+wSourceNames[t]);
        application.createModule<Quark>(wQuarkNames[t], quarkPar);
    }
    
    
    time_t t = time(0);
    struct tm* now = localtime(&t);
    
    char buffer[80];
    strftime(buffer, 80, "%Y%m%d",now);
    std::string current_date(buffer);
    strftime(buffer, 80, "%H%M%S",now);
    std::string current_time(buffer); 

    MContraction::TwoPion::Par twoPionPar;
    twoPionPar.output = "twopion/"+current_date+"/Full_" + current_time;
     
    twoPionPar.q0_noise_1     = "Qu1_0-1";
    twoPionPar.q0_noise_2   = "Qu1_0-2";
    twoPionPar.q_pos_noise_1     = "Qu1_p_1";
    twoPionPar.q_pos_noise_2    = "Qu1_p_2";
    twoPionPar.q_neg_noise_1     = "Qu1_n_1";
    twoPionPar.q_neg_noise_2     = "Qu1_n_2";
    twoPionPar.qs_pn    = "QS_PN";
    twoPionPar.qs_np    = "QS_NP";
    twoPionPar.v_qs_p0 = wQuarkNames;

    twoPionPar.mom    = positive_momentum;

    application.createModule<MContraction::TwoPion>("twoPion_Stochastic",
            twoPionPar);

    // execution
    application.saveParameterFile("twopion/"+current_date+"/xml/PiPi_"+current_time);
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
