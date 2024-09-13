
#include "EventAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern long int Nph[1152];
extern G4int Nevent;
extern std::ofstream PhEffect;

extern float Tph[1152][200000];
extern float Eph[1152][200000];
//extern G4double x0, yy0, z0, teta, phi, x, y, z, ux, uy, uz;
//extern G4double strx, stry, strz, strux, struy, struz;
//extern char fpartname[7];
extern char NazFile[1500];

long int q, p;

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(), fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{    
    //std::fill(Eph[0][0], Eph[1152][200000], 0);
    //std::fill(Tph[0][0], Tph[1152][200000], 0);
    std::fill_n(Nph, 1152, 0);

     for (p = 0; p < 1152; p++)
     {
         for(q = 0; q < 200000; q++)
         {
             Tph[p][q] = 0;
             Eph[p][q] = 0;
         }
     } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{   
    G4int NumSiPM, NumPhot;

    sprintf(NazFile, "D:\\MMH_reconst_data\\EVENTS_DATA\\EventM%05d.dat", Nevent);
    PhEffect.open(NazFile, std::ios_base::trunc);
    if(PhEffect.is_open())
    {
       PhEffect << "Energy" << "\t" << "Time" << "\t" << "Ncopy" << std::endl;
       for(NumSiPM = 0; NumSiPM < 1152; NumSiPM++)
       {
            for (NumPhot = 0; NumPhot < Nph[NumSiPM]; NumPhot++)
            {
                PhEffect << Eph[NumSiPM][NumPhot] << "\t" <<  Tph[NumSiPM][NumPhot] << "\t" << NumSiPM << std::endl;
            }
       }
        
    }
    //G4cout << "Number of event: " << Nevent << G4endl;
    PhEffect.close();
    Nevent++;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
