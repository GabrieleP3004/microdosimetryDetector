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

#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(AnalysisManager* pAnalysis)
{
	gun = new G4GeneralParticleSource();
	analysis = pAnalysis;
	numEv=0;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete gun;
}	

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	gun -> GeneratePrimaryVertex(anEvent);
	
	//print the number of events every NP events
	G4double NP=1000000.;
	numEv=numEv+1;
	G4double numD=numEv/NP;
	G4int numI=numEv/NP;
	if (numD==numI) {
		G4cout << "Number of events = " << numEv << G4endl;
	}

	//fill the firt ntuple
	G4double energy = gun -> GetParticleEnergy();
	analysis -> SetPrimaryEnergy(energy); 
}



