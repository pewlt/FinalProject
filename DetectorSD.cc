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
//
/// \file B4/B4c/src/CalorimeterSD.cc
/// \brief Implementation of the B4c::CalorimeterSD class

#include "DetectorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SensitiveDetector::SensitiveDetector(const G4String& name,
                             const G4String& hitsCollectionName,
                             G4int nofCells)
 : G4VSensitiveDetector(name),
   fNofCells(nofCells) //nombre de cellule = nombre de panel


{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection
    = new SlabHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

  // Create hits
  // fNofCells for cells + one more for total sums
  for (G4int i=0; i<fNofCells + 1; i++ ) {
//  for (G4int i=0; i<fNofCellsPannel; i++ ) {
    fHitsCollection->insert(new SlabHit());
      
  }
//  fHitsCollection->insert(new EventHit());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SensitiveDetector::ProcessHits(G4Step* step,
                                     G4TouchableHistory*)
{
  // energy deposit
  
  auto Edep = step->GetTotalEnergyDeposit();
  auto HitTime = step->GetPreStepPoint()->GetLocalTime();
  /*
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
  stepLength = step->GetStepLength();
  }

  if ( Edep==0. && stepLength == 0. ) return false;
*/
  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id
  //auto layerNumber = touchable->GetReplicaNumber(1);
  auto ChannelNbr = touchable->GetReplicaNumber();
  auto PannelNbr = touchable->GetReplicaNumber(1);

  // Get hit accounting data for this cell
  
auto hit = (*fHitsCollection)[PannelNbr];
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << PannelNbr;
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }
  
/*
 *  auto hit = (*fHitsCollection)[layerNumber]; //PannelNbr dans notre cas
  if ( ! hit ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << layerNumber;
    G4Exception("CalorimeterSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }
*/



  // Get hit for total accounting
//  auto hitTotal
//    = (*fHitsCollection)[fHitsCollection->entries()-1];

  // Add values
  hit->SetEdep(Edep);
  hit->SetHitTime(HitTime);
  hit->SetChannelNbr(ChannelNbr);
  hit->SetPanelNbr(PannelNbr);
  (*fHitsCollection)[4]->Add(Edep , HitTime);
//  hit->Add(Edep, HitTime);
//  hitTotal->Add(Edep, HitTime);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl
       << "-------->Hits Collection: in this event they are " << nofHits
       << " hits in the tracker chambers: " << G4endl;
     //for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
/*
  G4double fEdepTot = 0;
  G4double fTimeTot = 0;
  for (G4int i=0; i<fNofCells; i++ ) {
      auto hit = (*fHitsCollection)[i];
      fEdepTot += hit-> GetEdep();
      fTimeTot += hit-> GetHitTime();
      
  }
  fHitsCollection->insert(fEdepTot);
  fHitsCollection->insert(fTimeTot);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
