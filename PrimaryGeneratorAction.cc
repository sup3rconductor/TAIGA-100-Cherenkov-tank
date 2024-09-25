#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>

G4double theta, phi, ux, uy, uz, E0;
char fpartname[7];							
G4double ftheta, fphi, fEkin;
G4int fEvent, fpartnum;

G4double TankRad = 0.6 * m;
G4double z_mid, R_min, R_max, t_min, t_max;
G4double x_rand, y_rand, z_rand, R_rand, t_rand, x_up, y_up, z_up;
G4double X = 0. * m, Y = 0. * m, Z = 0. * m;
extern G4double Z0const;

extern std::ifstream rdata;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fParticleGun(0)
{
	G4int n_particle = 1;
	fParticleGun = new G4ParticleGun(n_particle);

	// default particle kinematic
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
	fParticleGun->SetParticleDefinition(particle);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
	fParticleGun->SetParticleEnergy(4. * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	/*//Circle radius where particle hits
	R_min = 0. * m;
	R_max = TankRad;
	R_rand = R_min + (R_max - R_min) * G4UniformRand();

	//Parameter for circle equation
	t_min = 0. * deg;
	t_max = 360. * deg;
	t_rand = t_min + (t_max - t_min) * G4UniformRand();

	//Coordinates for random directions
	x_rand = R_rand * cos(t_rand);
	y_rand = R_rand * sin(t_rand);
	z_rand = z_mid;

	//Scanning particle data from file
	rdata >> fEvent >> fpartname >> fpartnum >> ftheta >> fphi >> fEkin;

	std::cout << "Num of event: " << fEvent << "; " <<  "Particle name: " << fpartname << "; " << "Particle number: " << fpartnum << "; " << 
	"Theta angle: " << ftheta << "; " << "Phi angle: " << fphi << "; " << "E: " << fEkin << std::endl; 


	//Converting degrees to radians
	theta = ftheta * pi / 180.0;
	phi = fphi * pi / 180.0;

	//Setting type of particle to a particle gun
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle(fpartname);
	fParticleGun->SetParticleDefinition(particle); 

	//Launching position
	x_up = x_rand + 4500 * cm * sin(theta) * cos(phi);
	y_up = y_rand + 4500 * cm * sin(theta) * sin(phi);
	z_up = z_rand + 4500 * cm * cos(theta);

	fParticleGun->SetParticlePosition(G4ThreeVector(x_up, y_up, z_up));

	//Particle momentum direction
	ux = -sin(theta) * cos(phi);
	uy = -sin(theta) * sin(phi);
	uz = -cos(theta); 

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));*/

	//Default - vertical muon with 4 GeV
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
	fParticleGun->SetParticleDefinition(particle);

	fParticleGun->SetParticlePosition(G4ThreeVector(0.05 * m, 0.05 * m, 5. * m));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

	//Particle kinetic energy
	//fParticleGun->SetParticleEnergy(fEkin * GeV);
	fParticleGun->SetParticleEnergy(4. * GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


