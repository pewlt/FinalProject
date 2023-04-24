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
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorSD.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
    // Geometry parameters
  G4int fNofSlab = 16; //nombre de slab
  G4double SlabLength = 160.*cm; //longueur de la slab
  G4double SlabWidth = 10.*cm; //largeur de la slab
  G4double SlabThickness = 3.*cm; //epaisseur de la slab
  G4double gapThickness = 100.*cm; //ecart entre les deux plaques

  auto TotalThickness = 4 * SlabThickness + gapThickness; //epaisseur totale
  auto TotalWidth = fNofSlab * SlabWidth; //largeur totale
  auto worldSizeX = 1.2 * SlabLength; //longeur du monde primaire
  auto worldSizeY = 1.2 * TotalWidth; //largeur du monde primaire
  auto worldSizeZ = 1.2 * TotalThickness; //epaisseur du monde primaire
  
  //
  // Material
  //
   // Vacuum
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;

  new G4Material("Galactic", z=1., a=1.01*g/mole,density= CLHEP::universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Material* world_mat = nist->FindOrBuildMaterial("Galactic"); //matériel par défault du monde
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR"); //matériel par défaut de l'enveloppe
  G4Material* slab_mat = nist->FindOrBuildMaterial("G4_ANTHRACENE"); //matériel des slab
  /*
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
*/
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool fCheckOverlaps = true;
/*
  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking
*/




  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
//                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2);

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 world_mat,        // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps
  //
  // Envelope
  //
  auto EnvS = new G4Box("Envelope",                    // its name
    0.5 * SlabLength, 0.5 * TotalWidth, 2. * TotalThickness);  // its size

  auto EnvLV = new G4LogicalVolume(EnvS,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    EnvLV,                 // its logical volume
    "Envelope",               // its name
    worldLV,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking
/*
  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    shape1_mat,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  auto solidShape2 = new G4Trd("Shape2",  // its name
    0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
    0.5 * shape2_dz);  // its size

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    shape2_mat,                                        // its material
    "Shape2");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos2,                     // at position
    logicShape2,              // its logical volume
    "Shape2",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking
*/


  //
  //Slab
  //

  auto Slab
    = new G4Box("Slab",
                SlabLength/2, SlabWidth/2, SlabThickness/2);

//  auto SlabLV
    auto SlabLV
    = new G4LogicalVolume(
                 Slab,
//                 SLabS,           // its solid
                 slab_mat,  // its material
                 "SlabLV");         // its name

  //
  //Panel1
  //
  auto Panel1
    = new G4Box("Panel1",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel1LV
    = new G4LogicalVolume(
                 Panel1,
//                 layerS,           // its solid
                 slab_mat,  // its material
                 "Panel1");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0 , 0 , -(gapThickness/2 + SlabThickness/2 + SlabThickness)),          // at (0,0,0)
    Panel1LV,                  // its logical volume
    "Panel1",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel1",          // its name
                 SlabLV,          // its logical volume
                 Panel1LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel2
  auto Panel2
    = new G4Box("Panel2",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
    auto Panel2LV
    = new G4LogicalVolume(
                 Panel2,
//                 layerS,           // its solid
                 slab_mat,  // its material
                 "Panel2");         // its name

  
  G4RotationMatrix* r = new G4RotationMatrix;
  r ->rotateZ(M_PI/2.*rad);  
    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0 , 0 , -(gapThickness/2 + SlabThickness/2)),          // at (0,0,0)
    Panel2LV,                  // its logical volume
    "Panel2",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel2",          // its name
                 SlabLV,          // its logical volume
                 Panel2LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

  
  ////Panel3
  auto Panel3
    = new G4Box("Panel3",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel3LV
    = new G4LogicalVolume(
                 Panel3,
//                 layerS,           // its solid
                 slab_mat,  // its material
                 "Panel3");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0 , 0 , (gapThickness/2 + SlabThickness/2)),          // at (0,0,0)
    Panel3LV,                  // its logical volume
    "Panel3",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  


  new G4PVReplica(
                 "Panel3",          // its name
                 SlabLV,          // its logical volume
                 Panel3LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel4
  auto Panel4
    = new G4Box("Panel4",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
  auto Panel4LV
    = new G4LogicalVolume(
                 Panel4,
//                 layerS,           // its solid
                 slab_mat,  // its material
                 "Panel4");         // its name


    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0 , 0 , (gapThickness/2 + SlabThickness + SlabThickness/2)),          // at (0,0,0)
    Panel4LV,                  // its logical volume
    "Panel4",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel4",          // its name
                 SlabLV,          // its logical volume
                 Panel4LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica
  
  
 fScoringVolume = Panel4LV;
 
  // Set Shape2 as scoring volume
  //
 // fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  return worldPV;
}

void DetectorConstruction :: ConstructSD()
{
  G4int fNofPanel = 4; //nombre de panel
/*
  auto layerScintillatorSD
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD
    = new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);

*/    
  auto slabSD
    = new SensitiveDetector("SlabSD" , "SlabHitCollection" , fNofPanel); 
  G4SDManager::GetSDMpointer()->AddNewDetector(slabSD);
  SetSensitiveDetector("SlabLV" , slabSD);
    
    
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
