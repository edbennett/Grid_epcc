#include <Grid/Grid.h>

struct hmc_params {
  int starttraj;
  double beta;
  double mass;
  double hb_mu;
  double tlen;
  int nsteps;
};



hmc_params ReadCommandLineHMC(int argc, char** argv) {
  hmc_params HMCParams;
  if (Grid::GridCmdOptionExists(argv, argv + argc, "--starttraj")) {
     HMCParams.starttraj = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--starttraj"));
  } else {
    std::cout << Grid::GridLogError << "--starttraj must be specified" << std::endl;
  }

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

  if (Grid::GridCmdOptionExists(argv, argv + argc, "--hb_mu")) {
     HMCParams.hb_mu = std::stod(Grid::GridCmdOptionPayload(argv, argv + argc, "--hb_mu"));
  } else {
    std::cout << Grid::GridLogError << "--hb_mu must be specified" << std::endl;
  }

  return HMCParams;
}



int main(int argc, char **argv) {
  using namespace Grid;
  hmc_params HMCParams = ReadCommandLineHMC(argc, argv);
  
  typedef Representations<SpFundamentalRepresentation,
                          SpTwoIndexAntiSymmetricRepresentation>
      TheRepresentations;
  
  Grid_init(&argc, &argv);
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2>
      HMCWrapper;
  typedef SpWilsonTwoIndexAntiSymmetricImplR FermionImplPolicy;
  typedef SpWilsonTMTwoIndexAntiSymmetricFermionD FermionAction;
  typedef typename FermionAction::FermionField FermionField;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  HMCparameters HMCparams2;
  HMCparams2.StartTrajectory  = HMCParams.starttraj;
  HMCWrapper TheHMC(HMCparams2);

  TheHMC.Resources.AddFourDimGrid("gauge");

  // Checkpointer definition
  CheckpointerParameters CPparams;
  CPparams.config_prefix = "./hasenbusch/ckpoint_lat";
  CPparams.rng_prefix = "./hasenbusch/ckpoint_rng";
  CPparams.saveInterval = 1;
  CPparams.format = "IEEE64BIG";

  TheHMC.Resources.LoadNerscCheckpointer(CPparams);

  OneFlavourRationalParams Params(1.0e-6, 64.0, 2000, 1.0e-6, 16);

  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);

  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();

  RealD beta = HMCParams.beta;
  std::vector<RealD> hasenbusch({ HMCParams.hb_mu });
  SpWilsonGaugeActionR Waction(beta);


  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();

  SpTwoIndexAntiSymmetricRepresentation::LatticeField U(GridPtr);
 
  // LatticeGaugeField U(GridPtr);

  RealD mass = HMCParams.mass;
  Real pv_mass      = 10000.0;

  std::vector<Complex> boundary = {+1, +1, +1, -1};
  FermionAction::ImplParams bc(boundary);

  ConjugateGradient<FermionField> CG(1.0e-8, 2000, false);

  ActionLevel<HMCWrapper::Field, TheRepresentations> Level1(1);
  ActionLevel<HMCWrapper::Field, TheRepresentations> Level2(4);

  std::cout << GridLogMessage << "mass " << mass << std::endl;

  // Hasenbusch
  std::vector<RealD> light_den;
  std::vector<RealD> light_num;

  int n_hasenbusch = hasenbusch.size();
  light_den.push_back(0); // Target mass
  for(int h=0;h<n_hasenbusch;h++){
    light_den.push_back(hasenbusch[h]);
    light_num.push_back(hasenbusch[h]);
  }
  light_num.push_back(pv_mass);

  std::vector<FermionAction *> Numerators;
  std::vector<FermionAction *> Denominators;
  std::vector<TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy> *> Quotients;
  std::vector<OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy> *> Quotients2;

  // Untwisted, target mass
  Denominators.push_back(new FermionAction(U, *GridPtr, *GridRBPtr, mass, 0, bc));
  for(int h=0; h<n_hasenbusch; h++){
    std::cout << GridLogMessage << " 2f quotient Action  "<< light_num[h] << " / " << light_den[h]<< std::endl;
    Numerators.push_back  (new FermionAction(U, *GridPtr, *GridRBPtr, mass, hasenbusch[h], bc));
    Denominators.push_back(new FermionAction(U, *GridPtr, *GridRBPtr, mass, hasenbusch[h], bc));
  }
  // Quenched PV mass to balance the fractions
  Numerators.push_back(new FermionAction(U, *GridPtr, *GridRBPtr, pv_mass, 0, bc));

  for (int h=0; h < n_hasenbusch + 1; h++) {
    Quotients.push_back   (new TwoFlavourEvenOddRatioPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],CG,CG));
    Quotients2.push_back   (new OneFlavourEvenOddRatioRationalPseudoFermionAction<FermionImplPolicy>(*Numerators[h],*Denominators[h],Params));
  }

  for(int h=0;h<n_hasenbusch+1;h++){
    Level1.push_back(Quotients[h]);
    Level1.push_back(Quotients2[h]);
  }
  Level2.push_back(&Waction);
  
  TheHMC.TheAction.push_back(Level1);
  TheHMC.TheAction.push_back(Level2);

  TheHMC.Parameters.MD.MDsteps = HMCParams.nsteps;
  TheHMC.Parameters.MD.trajL = HMCParams.tlen;

  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();

  Grid_finalize();
}
