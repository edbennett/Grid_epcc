#include <Grid/Grid.h>

int main(int argc, char **argv) {
  using namespace Grid;
   
  typedef Representations< SpFundamentalRepresentation > TheRepresentations;

  Grid_init(&argc, &argv);
    
  typedef GenericSpHMCRunnerHirep<TheRepresentations, MinimumNorm2> HMCWrapper; // ok
  typedef SpWilsonImplR FermionImplPolicy;                                    // ok
  typedef SpWilsonFermionD FermionAction;                                     // ok??
  typedef SpWilsonTMFermionD HasenbuschAction;                                     // ok??
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
    
  RealD beta = 7.05;
    
  SpWilsonGaugeActionR Waction(beta);
    
  auto GridPtr = TheHMC.Resources.GetCartesian();
  auto GridRBPtr = TheHMC.Resources.GetRBCartesian();
    
  SpFundamentalRepresentation::LatticeField U(GridPtr);
    
  RealD mass = -0.81;
  RealD hasenbusch_twist = 0.3;

  FermionAction fermion(U, *GridPtr, *GridRBPtr, mass);
  HasenbuschAction numerator(U, *GridPtr, *GridRBPtr, mass, hasenbusch_twist);
  HasenbuschAction denominator(U, *GridPtr, *GridRBPtr, mass, hasenbusch_twist);

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
    
  TheHMC.Parameters.MD.MDsteps = 15;
  TheHMC.Parameters.MD.trajL   = 1.0;
    
  TheHMC.ReadCommandLine(argc, argv);
  TheHMC.Run();
    
  Grid_finalize();
}
