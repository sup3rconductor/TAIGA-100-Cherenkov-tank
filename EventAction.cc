
#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern long int Nph;
extern G4int Nevent;
extern std::ofstream PhEffect;
extern std::ofstream CHARGED;
char charged_file_path[1500];
char event_file_path[1500];

extern std::vector<double> Tph;
extern std::vector<double> Eph;
//extern G4double x0, yy0, z0, teta, phi, x, y, z, ux, uy, uz;
//extern G4double strx, stry, strz, strux, struy, struz;
//extern char fpartname[7];

long int q;

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(), fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{    
    Nph = 0;
    Eph.clear();
    Tph.clear();

    sprintf(charged_file_path, "E:\\TAIGA_SIMULATED_DATA\\CHARGED_MATRIX\\M03\\ChargedMatrix%05d.dat", Nevent);
    CHARGED.open(charged_file_path, std::ios_base::trunc);
    if (CHARGED.is_open()) 
    {
        CHARGED << "Particle" << "\t" << "Energy, MeV" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << std::endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{   
   G4int NumPhot;

    sprintf(event_file_path, "E:\\TAIGA_SIMULATED_DATA\\EVENTS_DATA\\M03\\Event%05d.dat", Nevent);
    PhEffect.open(event_file_path, std::ios_base::trunc);
    if(PhEffect.is_open())
    {
        PhEffect << "Energy" << "\t" << "Time" << std::endl;
        for(NumPhot = 0; NumPhot < Nph; NumPhot++)
        {
            PhEffect << Eph[NumPhot] << "\t" <<  Tph[NumPhot] << std::endl;
        }
        
    }
    G4cout << "Number of event: " << Nevent << G4endl;
    PhEffect.close();
    CHARGED.close();
    Nevent++; 
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
