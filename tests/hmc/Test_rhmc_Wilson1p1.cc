#include <string>
#include <Grid/Grid.h>


struct hmc_params {
  int save_freq;
  double beta;
  double m;
  double tlen;
  int nsteps;
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
  typedef HMCWrapperTemplate<PeriodicGimplD, MinimumNorm2, TheRepresentations> HMCWrapper; // Try to get double precision here too??
  // typedef GenericHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  typedef WilsonAdjImplD FermionImplPolicy;
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
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  //////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Collect actions, here use more encapsulation
  // need wrappers of the fermionic classes 
  // that have a complex construction
  // standard
  RealD beta = HMCParams.beta;
  WilsonGaugeActionD Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  // temporarily need a gauge field
  AdjointRepresentation::LatticeField U(GridPtr);

  Real mass = HMCParams.m;

  // Can we define an overloaded operator that does not need U and initialises
  // it with zeroes?
  // These lines are unecessary if BC are all periodic
  std::vector<Complex> boundary = {1,1,1,-1};
  FermionAction::ImplParams FParams(boundary);
  FermionAction FermOp(U, *GridPtr, *GridRBPtr, mass, FParams);

  // 1 flavour
  //OneFlavourRationalParams PfParams(1.0e-4, 64.0, 2000, 1.0e-10, 25, 100);
  OneFlavourRationalParams PfParams(1.0e-4, 64.0, 2000, 1.0e-6);
  OneFlavourEvenOddRationalPseudoFermionAction<FermionImplPolicy> WilsonNf1(FermOp, PfParams);

  //Smearing on/off
  WilsonNf1.is_smeared = false;

    // Collect actions
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  Level1.push_back(&WilsonNf1);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(5);
  Level2.push_back(&Waction);

  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);
  /////////////////////////////////////////////////////////////

  /*
    double rho = 0.1;  // smearing parameter
    int Nsmear = 2;    // number of smearing levels
    Smear_Stout<HMCWrapper::ImplPolicy> Stout(rho);
    SmearedConfiguration<HMCWrapper::ImplPolicy> SmearingPolicy(
        UGrid, Nsmear, Stout);
  */

  // HMC parameters are serialisable 
  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL   = HMCParams.tlen;

  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  TheHMC.Run();  // no smearing
  // TheHMC.Run(SmearingPolicy); // for smearing

  Grid_finalize();

} // main

