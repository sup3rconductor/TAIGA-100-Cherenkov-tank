#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Z0const, X0const, Y0const;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	/*	MATERIALS	*/

	G4double a, z;  //Atomic mass, atomic number
	G4double density, fractionmass;
	G4int ncomponents, nelements;


	//Chemical elements
	G4Element* elH = new G4Element("Hydrogen", "H", z = 1., a = 1.01 * g / mole);
	G4Element* elC = new G4Element("Carbon", "C", z = 6., a = 12.01 * g / mole);
	G4Element* elN = new G4Element("Nitrogen", "N", z = 7., a = 14.01 * g / mole);
	G4Element* elO = new G4Element("Oxygen", "O", z = 8., a = 16.00 * g / mole);
	G4Element* elSi = new G4Element("Silicium", "Si", z = 14., a = 28.09 * g / mole);
	G4Element* elAl = new G4Element("Aluminium", "Al", z = 13., a = 26.98 * g / mole);
	G4Element* elB = new G4Element("Boron", "B", z = 5., a = 10.812 * g / mole);
	G4Element* elFe = new G4Element("Ferrum", "Fe", z = 26., a = 55.85 * g / mole);
	G4Element* elF = new G4Element("Fluor", "F", z = 9., a = 18.99 * g / mole);
	G4Element* elNa = new G4Element("Sodium", "Na", z = 11., a = 22.99 * g / mole);
	G4Element* elHg = new G4Element("Mercury", "Hg", z = 80., a = 200.59 * g / mole);
	G4Element* elK = new G4Element("Kalium", "K", z = 19., a = 39.1 * g / mole);
	G4Element* elCa = new G4Element("Calcium", "Ca", z = 20., a = 40.08 * g / mole);

	//Air
	G4Material* Air = new G4Material("MAir", density = 1.290 * mg / cm3, ncomponents = 2);
	Air->AddElement(elN, fractionmass = 0.8);
	Air->AddElement(elO, fractionmass = 0.2);

	//Water
	G4Material* Water = new G4Material("Water", density = 1.0 * g/ cm3, ncomponents = 2);
	Water->AddElement(elH, nelements = 2);
	Water->AddElement(elO, nelements = 2);

	//Aluminium
	G4Material* AlMaterial = new G4Material("MAluminium", z = 13., a = 26.98 * g / mole, density = 2.7 * g / cm3);

	//Standart rock
	G4Material* Soil = new G4Material("Soil", z = 11., a = 22. * g / mole, density = 2.1 * g / cm3);

	//Concrete
	G4Material* Concrete = new G4Material("Concrete", density = 2.3 * g / cm3, ncomponents = 10);
  	Concrete->AddElement(elH, fractionmass = 0.01);
  	Concrete->AddElement(elO, fractionmass = 0.529);
  	Concrete->AddElement(elNa, fractionmass = 0.016);
  	Concrete->AddElement(elHg, fractionmass = 0.002);
  	Concrete->AddElement(elAl, fractionmass = 0.034);
  	Concrete->AddElement(elSi, fractionmass = 0.337);
  	Concrete->AddElement(elK, fractionmass = 0.013);
  	Concrete->AddElement(elCa, fractionmass = 0.044);
  	Concrete->AddElement(elFe, fractionmass = 0.014);
  	Concrete->AddElement(elC, fractionmass = 0.001);

	//Polyethylene
	G4Material* C2H4 = new G4Material("Polyethylene", density = 0.89 * g / cm3, ncomponents = 2);
	C2H4->AddElement(elC, nelements = 2);
	C2H4->AddElement(elH, nelements = 4);

	//BoronSilicate Glass
	G4Material* SiO2 = new G4Material("MSiO2", density = 2.1 * g / cm3, ncomponents = 2);
	SiO2->AddElement(elSi, nelements = 1);
	SiO2->AddElement(elO, nelements = 2);
	G4Material* B2O3 = new G4Material("MB2O3", density = 2.26 * g / cm3, ncomponents = 2);
	B2O3->AddElement(elB, nelements = 2);
	B2O3->AddElement(elO, nelements = 3);
	G4Material* Al2O3 = new G4Material("MAl2O3", density = 3.99 * g / cm3, ncomponents = 2);
	Al2O3->AddElement(elAl, nelements = 2);
	Al2O3->AddElement(elO, nelements = 3);
	G4Material* Na2O = new G4Material("MNa2O", density = 2.27 * g / cm3, ncomponents = 2);
	Na2O->AddElement(elAl, nelements = 2);
	Na2O->AddElement(elO, nelements = 1);

	G4Material* BoronSilicateGlass = new G4Material("PMTGlass", density = 2.5 * g / cm3, ncomponents = 4);
	BoronSilicateGlass->AddMaterial(SiO2, fractionmass = 80. * perCent);
	BoronSilicateGlass->AddMaterial(B2O3, fractionmass = 14. * perCent);
	BoronSilicateGlass->AddMaterial(Al2O3, fractionmass = 4. * perCent);
	BoronSilicateGlass->AddMaterial(Na2O, fractionmass = 2. * perCent);

	
	/*	OPTICAL PROPERTIES */


	/*	CHERENKOV WATER TANK	*/
	G4bool checkOverlaps = true;

	//World
	G4double world_sizeX = 7 * m;
	G4double world_sizeY = 7 * m;
	G4double world_sizeZ = 7 * m;

	G4Box* solidWorld = new G4Box("World_s", 0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");
	G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);


	//Soil mount dimensions
	G4double MountUpInnerRad = 0. * m, MountDownInnerRad = 0. * m;
	G4double MountUpOuterRad = 0.61 * m;
	G4double MountDownOuterRad = 2.61 * m;
	G4double MountHeight = 2. * m;
	G4double StartPhi = 0. * deg, StopPhi = 360. * deg;

	//Soil cylinder dimensions
	G4double CylInnerRad = 0. * m, CylOuterRad = 2.61 * m;
	G4double CylHeight = 1. * m;

	//Concrete tank dimensions
	G4double TankInnerRad = 0. * m, TankOuterRad = 0.610 * m;
	G4double TankHeight = CylHeight;
	G4double TankThickness = 0.01 * mm;

	//Polyethylene film dimensions
	G4double PEFilmInnerRad = 0. * m, PEFilmOuterRad = 0.6 * m;
	G4double PEFilmHeight = TankHeight - TankThickness;
	G4double PEFilmThickness = 0.2 * mm;

	//Water volume dimensions
	G4double WaterInnerRad = 0. * m, WaterOuterRad = PEFilmOuterRad - PEFilmThickness;
	G4double WaterHeight = PEFilmHeight - PEFilmThickness;

	//PMT glass dimensions 
	G4double UpperPartHAX = 133 * mm, UpperPartHAY = 133 * mm, UpperPartHAZ = 133 * mm;
	G4double UpperPartCut = 81.5 * mm, UpperPartDeposit = 44 * mm; 

	//Coordinates
	G4double XCyl = 0. * m, YCyl = 0. * m, ZCyl = 0. * m;												//Soil cylinder coordinates
	G4double XMount = XCyl, YMount = YCyl, ZMount = 0.5 * (CylHeight + MountHeight) + TankThickness;	//Soil mount coordinates
	G4double XTank = XCyl, YTank = YCyl, ZTank = ZCyl;													//Concrete tank coordinates
	G4double XPEFilm = XCyl, YPEFilm = YCyl, ZPEFilm = ZCyl;											//PE film coordinates
	G4double XWater = XCyl, YWater = YCyl, ZWater = ZCyl;												//Water volume coordinates
	G4double XUpperPartGlass = XCyl, YUpperPartGlass = YCyl, ZUpperPartGlass = 0.5 * WaterHeight - UpperPartHAZ + UpperPartDeposit;

	//Volumes
	G4Cons* solidSoilMount = {nullptr};

	G4Ellipsoid* solidUpperPartGlass = {nullptr};

	G4Tubs* solidSoilCylinder = {nullptr}, * solidTank = {nullptr}, * solidWaterVolume = {nullptr}, * solidPEFilm = {nullptr};

	G4LogicalVolume* logicSoilMount = {nullptr}, * logicSoilCylinder = {nullptr}, * logicTank = {nullptr}, * logicWaterVolume = {nullptr}, * logicPEFilm = {nullptr}, 
		* logicUpperPartGlass = {nullptr};

	G4VPhysicalVolume* physSoilMount = {nullptr}, * physSoilCylinder = {nullptr}, * physTank = {nullptr}, * physWaterVolume = {nullptr}, * physPEFilm = {nullptr}, 
		* physUpperPartGlass = {nullptr};

	//Building water tank
	solidSoilCylinder = new G4Tubs("Cyl_s", CylInnerRad, CylOuterRad, 0.5 * CylHeight, StartPhi, StopPhi);
	logicSoilCylinder = new G4LogicalVolume(solidSoilCylinder, Soil, "Cyl_l");
	physSoilCylinder = new G4PVPlacement(0, G4ThreeVector(XCyl, YCyl, ZCyl), logicSoilCylinder, "SOIL_CYLINDER", logicWorld, false, 0, checkOverlaps);

	solidSoilMount = new G4Cons("Mount_s", MountDownInnerRad, MountDownOuterRad, MountUpInnerRad, MountUpOuterRad, 0.5 * MountHeight, StartPhi, StopPhi);
	logicSoilMount = new G4LogicalVolume(solidSoilMount, Soil, "Cyl_l");
	physSoilMount = new G4PVPlacement(0, G4ThreeVector(XMount, YMount, ZMount), logicSoilMount, "SOIL_MOUNT", logicWorld, false, 0, checkOverlaps);

	solidTank = new G4Tubs("Tank_s", TankInnerRad, TankOuterRad, 0.5 * TankHeight, StartPhi, StopPhi);
	logicTank = new G4LogicalVolume(solidTank, Concrete, "Tank_l");
	physTank = new G4PVPlacement(0, G4ThreeVector(XTank, YTank, ZTank), logicTank, "CONCRETE_TANK", logicSoilCylinder, false, 0, checkOverlaps);

	solidPEFilm = new G4Tubs("PEFilm_s", PEFilmInnerRad, PEFilmOuterRad, 0.5 * PEFilmHeight, StartPhi, StopPhi);
	logicPEFilm = new G4LogicalVolume(solidPEFilm, C2H4, "PEFilm_l");
	physPEFilm = new G4PVPlacement(0, G4ThreeVector(XPEFilm, YPEFilm, ZPEFilm), logicPEFilm, "POLYETHYLENE_FILM", logicTank, false, 0, checkOverlaps);

	solidWaterVolume = new G4Tubs("Water_s", WaterInnerRad, WaterOuterRad, 0.5 * WaterHeight, StartPhi, StopPhi);
	logicWaterVolume = new G4LogicalVolume(solidWaterVolume, Water, "Water_l");
	physWaterVolume = new G4PVPlacement(0, G4ThreeVector(XWater, YWater, ZWater), logicWaterVolume, "WATER_VOLUME", logicPEFilm, false, 0, checkOverlaps);

	solidUpperPartGlass = new G4Ellipsoid("UpperPartGlass_s", UpperPartHAX, UpperPartHAY, UpperPartHAZ, -133 * mm, - UpperPartCut);
	logicUpperPartGlass = new G4LogicalVolume(solidUpperPartGlass, BoronSilicateGlass, "UpperPartGlass_l");
	physUpperPartGlass = new G4PVPlacement(0, G4ThreeVector(XUpperPartGlass, YUpperPartGlass, ZUpperPartGlass), logicUpperPartGlass, "UPPER_PART_GLASS", logicWaterVolume, false, 0, checkOverlaps);

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......