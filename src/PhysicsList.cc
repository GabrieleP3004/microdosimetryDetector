#include "PhysicsList.hh"

#include <G4EmLivermorePhysics.hh>
#include <G4DecayPhysics.hh> 
#include <G4ProductionCutsTable.hh>
#include <G4SystemOfUnits.hh>

#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsQGSP_BIC_HP.hh>
#include <G4IonBinaryCascadePhysics.hh>

#include <G4RegionStore.hh>

PhysicsList::PhysicsList()
{
	//Low energy EM physics 
	RegisterPhysics(new G4EmLivermorePhysics());
	
	//Decay physics
	RegisterPhysics(new G4DecayPhysics());
	
	//Hadronic physics
	RegisterPhysics(new G4HadronElasticPhysicsHP());
	RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
	RegisterPhysics(new G4IonBinaryCascadePhysics());
}


void PhysicsList::SetCuts()
{
	// default cut values
	defaultCutValue = 1.*micrometer;
	fCutForGamma     = defaultCutValue;
	fCutForElectron  = defaultCutValue;
	fCutForPositron  = defaultCutValue;
	fCutForProton    = defaultCutValue;
	
	// set cuts everywhere
	SetCutValue(fCutForGamma, "gamma");
	SetCutValue(fCutForElectron, "e-");
	SetCutValue(fCutForPositron, "e+");
	SetCutValue(fCutForProton, "proton");
  
	
	// set different cuts in the high precision region
	G4String regName = "highPRegion";
	G4Region* region = G4RegionStore::GetInstance()->GetRegion(regName);
	
	G4ProductionCuts* cuts = new G4ProductionCuts;
	G4double highPCut = 0.1*um;
	
	cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("gamma"));
	cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("e-"));
	cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("e+"));
	cuts->SetProductionCut(highPCut, G4ProductionCuts::GetIndex("proton"));
	
	region->SetProductionCuts(cuts);
	
	// Update the production cuts table energy range
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250.*eV,10.*GeV);  
}
