//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

AnalysisManager::AnalysisManager() 
{
	factoryOn = false;

	// Initialization ntuple
	for (G4int k=0; k<MaxNtCol; k++) fNtColId[k] = 0;
}

AnalysisManager::~AnalysisManager() 
{
}

void AnalysisManager::book() 
{ 
	G4AnalysisManager* manager = G4AnalysisManager::Instance();
  
	manager->SetVerboseLevel(1);
 
	// Create a root file
	G4String fileName = "output.csv";

	// Create directories (not supported by csv)  
	//manager->SetNtupleDirectoryName("output Folder");
  

	G4bool fileOpen = manager->OpenFile(fileName);
	if (!fileOpen) {
		G4cout << "\n---> HistoManager::book(): cannot open " 
			<< fileName << G4endl;
		return;
	}

	manager->SetFirstNtupleId(1);

/*	//Create Primary Energy Ntuple
	manager -> CreateNtuple("1_primary_energy", "primary_energy");
	fNtColId[0] = manager -> CreateNtupleDColumn("E_keV");
	manager -> FinishNtuple();
*/
	// Create 4 series of Ntuple for the 4 detectors
	G4double SVside[4] = { 50, 300., 100., 200.};
	for (int i=0; i<4; i++)
	{
		
		//Create Energy Deposition within SV Ntuple
		std::ostringstream ntName1; ntName1 << SVside[i] << "um_deposited_energy";
		manager -> CreateNtuple(ntName1.str(), "deposited_energy");
		fNtColId[0+5*i] = manager -> CreateNtupleDColumn("E_keV");	//deposited energy 
		fNtColId[1+5*i] = manager -> CreateNtupleDColumn("l_um");	//chord length
 
		manager -> FinishNtuple();

		// creating Energy lost by primaries Ntuple
		std::ostringstream ntName2; ntName2 << SVside[i] << "um_energy_released";
		manager -> CreateNtuple(ntName2.str(), "energy_released");
		fNtColId[2+5*i] = manager -> CreateNtupleDColumn("Elost_keV");
		fNtColId[3+5*i] = manager -> CreateNtupleDColumn("Ein_keV");
		fNtColId[4+5*i] = manager -> CreateNtupleDColumn("Eout_keV");
		manager -> FinishNtuple();
	}

/*	//creating a ntuple, containing the information about secondary particles
	manager -> CreateNtuple("3_secondary_particles", "secondary_particles");
	fNtColId[3] = manager -> CreateNtupleDColumn("AA");
	fNtColId[4] = manager -> CreateNtupleDColumn("ZZ");
	fNtColId[5] = manager -> CreateNtupleDColumn("Kin_keV");
	manager -> FinishNtuple();
*/
  
	factoryOn = true;    
}


/*void AnalysisManager::SetPrimaryEnergy(G4double energy)
{
	G4AnalysisManager* manager = G4AnalysisManager::Instance();
	manager -> FillNtupleDColumn(1, fNtColId[0], energy/keV);
	manager -> AddNtupleRow(1); 
}
*/

void AnalysisManager::StoreEnergyDeposition(G4double edep, G4double pathlen, G4int SVid)
{
	G4AnalysisManager* manager = G4AnalysisManager::Instance();
	manager -> FillNtupleDColumn(1+2*SVid, fNtColId[0+5*SVid], edep/keV);
	manager -> FillNtupleDColumn(1+2*SVid, fNtColId[1+5*SVid], pathlen/um);
	manager -> AddNtupleRow(1+2*SVid); 
}

void AnalysisManager::StorePrimaryEnergyLost(G4double elost, G4double ein, G4double eout, G4int SVid)
{
	G4AnalysisManager* manager = G4AnalysisManager::Instance();
	manager -> FillNtupleDColumn(2+2*SVid, fNtColId[2+5*SVid], elost/keV);
	manager -> FillNtupleDColumn(2+2*SVid, fNtColId[3+5*SVid], ein/keV);
	manager -> FillNtupleDColumn(2+2*SVid, fNtColId[4+5*SVid], eout/keV);
	manager -> AddNtupleRow(2+2*SVid);
}

/*void AnalysisManager::FillSecondaries(G4int AA, G4double charge, G4double energy)
  {

    G4AnalysisManager* manager = G4AnalysisManager::Instance();
    manager -> FillNtupleDColumn(3, fNtColId[3], AA);
    manager -> FillNtupleDColumn(3, fNtColId[4], charge);
    manager -> FillNtupleDColumn(3, fNtColId[5], energy/keV);
    manager -> AddNtupleRow(3);  
  }
*/
 
void AnalysisManager::finish() 
{   
	if (factoryOn) 
	{
		G4AnalysisManager* manager = G4AnalysisManager::Instance();    
		manager -> Write();
		manager -> CloseFile();  
      
		delete G4AnalysisManager::Instance();
		factoryOn = false;
	}
}












