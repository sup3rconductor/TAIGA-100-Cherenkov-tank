 
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern G4int Nevent;
extern G4double Z0const;
extern char fpartname[7];
G4int copyNo, strCopy;
long int Nph[1152];
float Tph[1152][200000];
float Energy, Epart, Eph[1152][200000];
char NazFileMU[50], NazFileChrgd[50];
extern std::ofstream MUON, CHARGED;
//G4double strx, stry, strz, strux, struy, struz;

SteppingAction::SteppingAction(EventAction* eventAction)
    : G4UserSteppingAction(), fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* theTrack = step->GetTrack();
  G4String Name = theTrack->GetDefinition()->GetParticleName();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4StepPoint* prePoint = step->GetPreStepPoint();
  G4String vname = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName(); 



  if (vname == "sipm_l") // SiPM volume
  {
    if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      Energy = theTrack->GetKineticEnergy()/eV;
      copyNo = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();
      if(Energy>=1.9 && Energy<=4.0 && Nph[copyNo]< 200000)
      {
        Eph[copyNo][Nph[copyNo]] = Energy;                             //Energy of Num N photon in Num copyNo of PMT
        Tph[copyNo][Nph[copyNo]] = theTrack->GetGlobalTime()/ns;        //Time of Num N photon in Num copyNo of PMT
        Nph[copyNo]++;                                                  //Plus one photon registered by Num copyNo of PMT
      }
    // collect optical photons 
      theTrack->SetTrackStatus(fStopAndKill);
    }   
  }

  if (vname == "strip_l" || vname == "tdlr_l" || vname == "exPS_l" || vname == "shell_l" || vname == "hollow_l" || vname == "glue_l")
  {
      strCopy = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();

      if (Name == "mu+" || Name == "mu-")
      {
          //Epart = theTrack()->GetKineticEnergy()/MeV;
          if (MUON.is_open()) MUON << Nevent << "\t" << Name.c_str() << "\t" << strCopy << "\t" << theTrack->GetKineticEnergy() / MeV << "\t" << prePoint->GetPosition().x() << "\t" <<
              prePoint->GetPosition().y() << "\t" << prePoint->GetPosition().z() << std::endl;
      }

      if (Name == "mu+" || Name == "mu-" || Name == "e+" || Name == "e-" || Name == "proton" || Name == "pi+" || Name == "pi-")
      {
          //Epart = theTrack()->GetKineticEnergy() / MeV;
          if (CHARGED.is_open()) CHARGED << Nevent << "\t" << Name.c_str() << "\t" << strCopy << "\t" << theTrack->GetKineticEnergy() / MeV << "\t" << prePoint->GetPosition().x() << "\t" <<
              prePoint->GetPosition().y() << "\t" << prePoint->GetPosition().z() << std::endl;
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

