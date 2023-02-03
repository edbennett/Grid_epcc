/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: 

Copyright (C) 2015-2016

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: Guido Cossu
Author: David Murphy

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
#include <Grid/Grid.h>

#ifdef GRID_DEFAULT_PRECISION_DOUBLE
#define MIXED_PRECISION
#endif

struct hmc_params {
  int save_freq;
  double beta;
  double m;
  double tlen;
  int nsteps;
  std::string serial_seed = "1 2 3 4 5";
  std::string parallel_seed = "6 7 8 9 10";
  int Ls = 8;
  double pv_mass = 1.8;
  double M5 = 1.8; // domain wall height
  double b = 1.0; // controls exactly what action is being used
  double c = 0.0; // b-c=1.0 (fixed?); alpha=b+c=1 is Shamir (see e.g. 1411.7017)
  std::string cnfg_dir = ".";
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
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--serialseed")) {
    HMCParams.serial_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--serialseed");
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--parallelseed")) {
    HMCParams.parallel_seed = Grid::GridCmdOptionPayload(argv, argv + argc, "--parallelseed");
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--pv_mass")) {
    HMCParams.pv_mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--pv_mass"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--dw_mass")) {
    HMCParams.M5 = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--dw_mass"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mobius_b")) {
    HMCParams.b = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mobius_b"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mobius_c")) {
    HMCParams.c = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mobius_c"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--Ls")) {
    HMCParams.Ls = std::stoi(Grid::GridCmdOptionPayload(argv, argv + argc, "--Ls"));
  }
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--cnfg_dir")) {
    HMCParams.cnfg_dir = Grid::GridCmdOptionPayload(argv, argv + argc, "--cnfg_dir");
  }

  return HMCParams;
}


int main(int argc, char **argv) {
  using namespace Grid;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);

  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  // here make a routine to print all the relevant information on the run
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

   // Typedefs to simplify notation
  typedef Representations<FundamentalRepresentation, AdjointRepresentation> TheRepresentations;
  typedef WilsonAdjImplR FermionImplPolicy;
  typedef MobiusFermion<WilsonAdjImplR> FermionAction;
  typedef MobiusEOFAFermion<WilsonAdjImplR> FermionEOFAAction;
  typedef typename FermionAction::FermionField FermionField;

  typedef Grid::XmlReader       Serialiser;
  
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //  typedef GenericHMCRunner<LeapFrog> HMCWrapper;
  //  MD.name    = std::string("Leap Frog");
  typedef GenericHMCRunnerHirep<TheRepresentations, ForceGradient> HMCWrapper;

  HMCWrapper TheHMC;
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL = HMCParams.tlen;

  // Grid from the command line arguments --grid and --mpi
  TheHMC.Resources.AddFourDimGrid("gauge"); // use default simd lanes decomposition
  
  CheckpointerParameters CPparams;
  CPparams.config_prefix = HMCParams.cnfg_dir + "/ckpoint_EODWF_lat";
  CPparams.rng_prefix    = HMCParams.cnfg_dir + "/ckpoint_EODWF_rng";
  CPparams.saveInterval  = HMCParams.save_freq;
  CPparams.format        = "IEEE64BIG";
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

  const int Ls      = HMCParams.Ls;
  Real beta         = HMCParams.beta;
  Real oneflavour_mass = HMCParams.m;
  Real pv_mass      = HMCParams.pv_mass;
  RealD M5  = HMCParams.M5;
  RealD b   = HMCParams.b;
  RealD c   = HMCParams.c;

  auto GridPtr   = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
  auto FGrid     = SpaceTimeGrid::makeFiveDimGrid(Ls,GridPtr);
  auto FrbGrid   = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,GridPtr);

  WilsonGaugeActionR GaugeAction(beta);

  // temporarily need a gauge field
  AdjointRepresentation::LatticeField U(GridPtr);

  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams FParams(boundary);

  double ActionStoppingCondition     = 1e-10;
  double DerivativeStoppingCondition = 1e-6;
  double MaxCGIterations = 30000;

  ////////////////////////////////////
  // Collect actions
  ////////////////////////////////////
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(8);

  ////////////////////////////////////
  // OneFlavour action
  ////////////////////////////////////

  // DJM: setup for EOFA ratio (Mobius)
  OneFlavourRationalParams OFRp;
  OFRp.lo       = 1.0e-4;
  OFRp.hi       = 64.0;
  OFRp.MaxIter  = 10000;
  OFRp.tolerance= 1.0e-9;
  OFRp.degree   = 20;
  OFRp.precision= 100;

  
  FermionEOFAAction OneFlavour_Op_L (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , oneflavour_mass, oneflavour_mass, pv_mass, 0.0, -1, M5, b, c, FParams);
  FermionEOFAAction OneFlavour_Op_R (U , *FGrid , *FrbGrid , *GridPtr , *GridRBPtr , pv_mass, oneflavour_mass,      pv_mass, -1.0, 1, M5, b, c, FParams);

  ConjugateGradient<FermionField>      ActionCG(ActionStoppingCondition,MaxCGIterations);
  ConjugateGradient<FermionField>  DerivativeCG(DerivativeStoppingCondition,MaxCGIterations);
  ExactOneFlavourRatioPseudoFermionAction<FermionImplPolicy>
    EOFA(OneFlavour_Op_L, OneFlavour_Op_R,
	 ActionCG,
	 ActionCG, ActionCG,
	 DerivativeCG, DerivativeCG, 
	 OFRp, true);
  Level1.push_back(&EOFA);


  /////////////////////////////////////////////////////////////
  // Gauge action
  /////////////////////////////////////////////////////////////
  Level2.push_back(&GaugeAction);
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.ReadCommandLine(argc, argv);
  std::cout << GridLogMessage << " Action complete "<< std::endl;

  /////////////////////////////////////////////////////////////
  // HMC parameters are serialisable

  std::cout << GridLogMessage << " Running the HMC "<< std::endl;
  TheHMC.Run();  // no smearing

  Grid_finalize();
} // main



