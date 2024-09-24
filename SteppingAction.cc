 
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

extern G4int Nevent;
extern G4double Z0const;
extern char fpartname[7];
G4int copyNo, strCopy;
long int Nph;
std::vector<double> Tph, Eph;
float Energy, Epart;
extern std::ofstream CHARGED;
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



  if (vname == "UpperPartPhotCath_l" || vname == "DownPartPhotCath_l") //Photocathode volume
  {
    if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      Energy = theTrack->GetKineticEnergy()/eV;
      copyNo = prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo();
      if(Energy >= 1.9 && Energy <= 3.3 && Nph< 200000)
      {
        Eph.push_back(Energy);                             //Energy of Num N photon in Num copyNo of PMT
        Tph.push_back(theTrack->GetGlobalTime()/ns);        //Time of Num N photon in Num copyNo of PMT
        Nph++;                                                  //Plus one photon registered by Num copyNo of PMT
      }
    // collect optical photons 
      theTrack->SetTrackStatus(fStopAndKill);
    }   
  }

  if (vname == "Tank_l" || vname == "PEFilm_l" || vname == "Water_l")
  {
      if (Name == "mu+" || Name == "mu-" || Name == "e+" || Name == "e-" || Name == "proton" || Name == "pi+" || Name == "pi-")
      {
          //Epart = theTrack()->GetKineticEnergy() / MeV;
          if (CHARGED.is_open()) CHARGED << Name.c_str() << "\t" << theTrack->GetKineticEnergy() / MeV << "\t" << prePoint->GetPosition().x() << "\t" <<
              prePoint->GetPosition().y() << "\t" << prePoint->GetPosition().z() << std::endl;
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

