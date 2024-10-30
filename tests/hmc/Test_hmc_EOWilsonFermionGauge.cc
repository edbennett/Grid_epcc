/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_hmc_WilsonFermionGauge.cc

Copyright (C) 2016

Author: Guido Cossu <guido.cossu@ed.ac.uk>

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
directory
*************************************************************************************/
/*  END LEGAL */
#include <string>
#include <Grid/Grid.h>

struct hmc_params {
  int save_freq;
  double beta;
  double m;
  double tlen;
  int nsteps;
  int nPV;
  double mPV;
  std::string serial_seed = "1 2 3 4 5";
  std::string parallel_seed = "6 7 8 9 10";
};

hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--savefreq")) {
    HMCParams.save_freq = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--savefreq"));
  } else {
    std::cout << Grid::GridLogError << "--savefreq must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--beta")) {
    HMCParams.beta = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--beta"));
  } else {
    std::cout << Grid::GridLogError << "--beta must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--fermionmass")) {
    HMCParams.m = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--fermionmass"));
  } else {
    std::cout << Grid::GridLogError << "--fermionmass must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--tlen")) {
    HMCParams.tlen = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--tlen"));
  } else {
    std::cout << Grid::GridLogError << "--tlen must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--nsteps")) {
    HMCParams.nsteps = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--nsteps"));
  } else {
    std::cout << Grid::GridLogError << "--nsteps must be specified" << std::endl;
    exit(1);
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--npv")) {
    HMCParams.nPV = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--npv"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mpv")) {
    HMCParams.mPV = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mpv"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--serialseed")) {
    HMCParams.serial_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--serialseed");
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--parallelseed")) {
    HMCParams.parallel_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--parallelseed");
  }
  return HMCParams;
}


int main(int argc, char **argv) {
  using namespace Grid;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);

  // Here change the allowed (higher) representations
  typedef Representations< FundamentalRepresentation, AdjointRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonAdjImplR FermionImplPolicy;
  typedef WilsonAdjFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  HMCWrapper TheHMC;

  // Grid from the command line
  TheHMC.Resources.AddFourDimGrid("gauge");
  // Possibile to create the module by hand 
  // hardcoding parameters or using a Reader


  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "cnfg/ckpoint_lat";
  CPparams.rng_prefix = "rand/ckpoint_rng";
  CPparams.saveInterval = HMCParams.save_freq;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = HMCParams.serial_seed;
  RNGpar.parallel_seeds = HMCParams.parallel_seed;
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  // here there is too much indirection 
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  RealD beta = HMCParams.beta;
  WilsonGaugeActionR Waction(beta);
  
  // temporarily need a gauge field
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  AdjointRepresentation::LatticeField U(GridPtr);
  LatticeGaugeField Ufun(GridPtr);

  Real mass = HMCParams.m;

  // Can we define an overloaded operator that does not need U and initialises
  // it with zeroes?
  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  TwoFlavourEvenOddPseudoFermionAction<FermionImplPolicy> Nf2(FermOp, CG, CG);

    // Set smearing (true/false), default: false
  Nf2.is_smeared = false;


  int N_PV = HMCParams.nPV;
  ConjugateGradient<WilsonFermionD::FermionField> CG_PV(1.0e-8, 2000, false);
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplD>> PVs;
  for (int PV_idx = 0; PV_idx < N_PV; PV_idx++ ) {
    double quenchedMass = 10.0;
    double PVMass = HMCParams.mPV;
    WilsonFermionD NumOp(Ufun, *GridPtr, *GridRBPtr, PVMass);
    WilsonFermionD DenOp(Ufun, *GridPtr, *GridRBPtr, quenchedMass);
    TwoFlavourEvenOddRatioPseudoFermionAction<WilsonImplD> PV(NumOp, DenOp, CG_PV, CG_PV);
    PV.is_smeared = true;
    PVs.push_back(PV);
  }

  // Collect actions
  ActionLevel<LatticeGaugeField, TheRepresentations> Level1(1);
  Level1.push_back(&Nf2);

  ActionLevel<LatticeGaugeField, TheRepresentations> Level2(4);
  Level2.push_back(&Waction);

  for (int PV_idx = 0; PV_idx < N_PV / 4; PV_idx++) {
    Level1.push_back(&(PVs[PV_idx]));
  }

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////


  double rho = 0.1;  // smearing parameter
  int Nsmear = 3;    // number of smearing levels
  Smear_Stout<HMCWrapper::ImplPolicy> Stout(rho);
  SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(
        GridPtr, Nsmear, Stout);

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL   = HMCParams.tlen;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main
