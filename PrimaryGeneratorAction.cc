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

G4double theta, phi, ux, uy, uz, E0;
char fpartname[7];							//���������������, ���� ������ � ������ �� ���������, ��������� ������ ��������� �� ������������ �����
G4int fEvent, fpartnum;
G4double ftheta, fphi, fEkin;


G4double Length = 2. * m, Width = 2. * m;
G4double x_min, x_max, y_min, y_max, z_mid;
G4double x_rand, y_rand, z_rand, x_up, y_up, z_up;
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
	//Maximal and minimal X coordinate
	x_min = X - 0.5 * Length;
	x_max = X + 0.5 * Length;

	//Maximal and minimal Y coordinate
	y_min = Y - 0.5 * Width;
	y_max = Y + 0.5 * Width;

	//Z coordinate of center
	z_mid = Z;

	//Coordinates for random directions
	x_rand = x_min + (x_max - x_min) * G4UniformRand();
	y_rand = y_min + (y_max - y_min) * G4UniformRand();
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
	x_up = x_rand + 500 * cm * sin(theta) * cos(phi);
	y_up = y_rand + 500 * cm * sin(theta) * sin(phi);
	z_up = z_rand + 500 * cm * cos(theta);

	fParticleGun->SetParticlePosition(G4ThreeVector(x_up, y_up, z_up));

	//Particle momentum direction
	ux = -sin(theta) * cos(phi);
	uy = -sin(theta) * sin(phi);
	uz = -cos(theta); 

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));

	//Particle kinetic energy
	fParticleGun->SetParticleEnergy(fEkin * GeV);
	fParticleGun->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


