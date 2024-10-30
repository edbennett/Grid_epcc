#include <Grid/Grid.h>

struct hmc_params {
  double beta;
  double mass;
  double hb_mass;
  double tlen;
  int nsteps;
};



hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--nsteps")) {
     HMCParams.nsteps = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--nsteps"));
  } else {
    std::cout << Grid::GridLogError << "--nsteps must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--tlen")) {
     HMCParams.tlen = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--tlen"));
  } else {
    std::cout << Grid::GridLogError << "--tlen must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--beta")) {
     HMCParams.beta = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--beta"));
  } else {
    std::cout << Grid::GridLogError << "--beta must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--mass")) {
     HMCParams.mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--mass"));
  } else {
    std::cout << Grid::GridLogError << "--mass must be specified" << std::endl;
  }

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--hb_mass")) {
     HMCParams.hb_mass = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--hb_mass"));
  } else {
    std::cout << Grid::GridLogError << "--hb_mass must be specified" << std::endl;
  }

  return HMCParams;
}


int main(int argc, char **argv) {
  using namespace Grid;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);
   
  typedef Representations< SpFundamentalRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
    
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper; // ok
  typedef SpWilsonImplR FermionImplPolicy;                                    // ok
  typedef SpWilsonTMFermionD FermionAction;                                     // ok??
  typedef typename FermionAction::FermionField FermionField;                  // ok?
 
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   
  HMCWrapper TheHMC;
    
  TheHMC.Resources.AddFourDimGrid("gauge");
    
  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";
    
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
    
  RealD beta = HMCParams.beta;
    
  SpWilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
    
  SpFundamentalRepresentation::LatticeField U(GridPtr);
    
  RealD mass = HMCParams.mass;
  RealD hasenbusch_twist = HMCParams.hb_mass;

  std::vector<Complex> boundary = {+1, +1, +1, -1};
  FermionAction::ImplParams bc(boundary);

  FermionAction fermion(U, *GridPtr, *GridRBPtr, mass, 0, bc);
  FermionAction numerator(U, *GridPtr, *GridRBPtr, mass, hasenbusch_twist, bc);
  FermionAction denominator(U, *GridPtr, *GridRBPtr, mass, hasenbusch_twist, bc);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  TwoFlavourRatioPseudoFermionAction<FermionImplPolicy> quotient(numerator, fermion, CG, CG);
  TwoFlavourPseudoFermionAction<FermionImplPolicy> remainder(denominator, CG, CG);

  quotient.is_smeared = false;
  remainder.is_smeared = false;
  
  ActionLevel<HMCWrapper::Field, TheRepresentations > Level1(1);
  Level1.push_back(&quotient);
    
  ActionLevel<HMCWrapper::Field, TheRepresentations > Level2(2);
  Level2.push_back(&remainder);

  ActionLevel<HMCWrapper::Field, TheRepresentations > Level3(2);
  Level3.push_back(&Waction);
    
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  TheHMC.TheAction.push_back(Level3);
    
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL   = HMCParams.tlen;
    
  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();
    
  Grid_finalize();
}
