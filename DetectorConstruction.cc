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
#include "G4NistManager.hh"

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

	//Vacuum inside PMT
	G4NistManager * NIST_man = G4NistManager::Instance();
	G4Material * Vacuum = NIST_man->FindOrBuildMaterial("G4_Galactic");


	
	/*	OPTICAL PROPERTIES */

	//Distillated water optical properties
	std::vector<G4double> PhotonEnergyWater = {
	  2.034 * eV, 2.068 * eV, 2.103 * eV, 2.139 * eV, 2.177 * eV, 2.216 * eV,
	  2.256 * eV, 2.298 * eV, 2.341 * eV, 2.386 * eV, 2.433 * eV, 2.481 * eV,
	  2.532 * eV, 2.585 * eV, 2.640 * eV, 2.697 * eV, 2.757 * eV, 2.820 * eV,
	  2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV,
	  3.353 * eV, 3.446 * eV, 3.545 * eV, 3.649 * eV, 3.760 * eV, 3.877 * eV,
	  4.002 * eV, 4.136 * eV };

	std::vector<G4double> RefractiveIndexWater = {
	  1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
	  1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
	  1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356,
	  1.3568, 1.3572, 1.358,  1.3585, 1.359,  1.3595, 1.36,   1.3608 };

	std::vector<G4double> AbsorptionLengthWater = {
	  3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
	  15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
	  45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
	  52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m,
	  30.000 * m, 28.500 * m, 27.000 * m, 24.500 * m, 22.000 * m, 19.500 * m,
	  17.500 * m, 14.500 * m };

	/*std::vector<G4double> energy_water = {
	  1.56962 * eV, 1.58974 * eV, 1.61039 * eV, 1.63157 * eV, 1.65333 * eV,
	  1.67567 * eV, 1.69863 * eV, 1.72222 * eV, 1.74647 * eV, 1.77142 * eV,
	  1.7971 * eV,  1.82352 * eV, 1.85074 * eV, 1.87878 * eV, 1.90769 * eV,
	  1.93749 * eV, 1.96825 * eV, 1.99999 * eV, 2.03278 * eV, 2.06666 * eV,
	  2.10169 * eV, 2.13793 * eV, 2.17543 * eV, 2.21428 * eV, 2.25454 * eV,
	  2.29629 * eV, 2.33962 * eV, 2.38461 * eV, 2.43137 * eV, 2.47999 * eV,
	  2.53061 * eV, 2.58333 * eV, 2.63829 * eV, 2.69565 * eV, 2.75555 * eV,
	  2.81817 * eV, 2.88371 * eV, 2.95237 * eV, 3.02438 * eV, 3.09999 * eV,
	  3.17948 * eV, 3.26315 * eV, 3.35134 * eV, 3.44444 * eV, 3.54285 * eV,
	  3.64705 * eV, 3.75757 * eV, 3.87499 * eV, 3.99999 * eV, 4.13332 * eV,
	  4.27585 * eV, 4.42856 * eV, 4.59258 * eV, 4.76922 * eV, 4.95999 * eV,
	  5.16665 * eV, 5.39129 * eV, 5.63635 * eV, 5.90475 * eV, 6.19998 * eV };

	std::vector<G4double> mie_water = {
		167024.4 * m, 158726.7 * m, 150742 * m,   143062.5 * m, 135680.2 * m,
		128587.4 * m, 121776.3 * m, 115239.5 * m, 108969.5 * m, 102958.8 * m,
		97200.35 * m, 91686.86 * m, 86411.33 * m, 81366.79 * m, 76546.42 * m,
		71943.46 * m, 67551.29 * m, 63363.36 * m, 59373.25 * m, 55574.61 * m,
		51961.24 * m, 48527.00 * m, 45265.87 * m, 42171.94 * m, 39239.39 * m,
		36462.50 * m, 33835.68 * m, 31353.41 * m, 29010.30 * m, 26801.03 * m,
		24720.42 * m, 22763.36 * m, 20924.88 * m, 19200.07 * m, 17584.16 * m,
		16072.45 * m, 14660.38 * m, 13343.46 * m, 12117.33 * m, 10977.70 * m,
		9920.416 * m, 8941.407 * m, 8036.711 * m, 7202.470 * m, 6434.927 * m,
		5730.429 * m, 5085.425 * m, 4496.467 * m, 3960.210 * m, 3473.413 * m,
		3032.937 * m, 2635.746 * m, 2278.907 * m, 1959.588 * m, 1675.064 * m,
		1422.710 * m, 1200.004 * m, 1004.528 * m, 833.9666 * m, 686.1063 * m };


	// ### OpNovice
	G4double cosTeta = 0.99;      // c
	G4double cosTetaBack = 0.99;  // c
	G4double f_ratio = 0.8;       // ###*/

	G4MaterialPropertiesTable* WaterOptProperties = new G4MaterialPropertiesTable();
	WaterOptProperties->AddProperty("RINDEX", PhotonEnergyWater, RefractiveIndexWater, PhotonEnergyWater.size());
	WaterOptProperties->AddProperty("ABSLENGTH", PhotonEnergyWater, AbsorptionLengthWater, PhotonEnergyWater.size());

	/*WaterOptProperties->AddProperty("MIEHG", energy_water, mie_water, mie_water.size());  // ### OpNovice
	WaterOptProperties->AddConstProperty("MIEHG_FORWARD", cosTeta);
	WaterOptProperties->AddConstProperty("MIEHG_BACKWARD", cosTeta);
	WaterOptProperties->AddConstProperty("MIEHG_FORWARD_RATIO", f_ratio);*/

	Water->SetMaterialPropertiesTable(WaterOptProperties);

	//BoronSilicate Glass optical properties
	G4double PhotonEnergy[2] = {1.9 * eV, 4.0 * eV};
	G4double RefractiveIndexGlass[2] = {1.5, 1.5};
	G4double AbsorptionLengthGlass[2] = {40. * m, 40. * m};

	G4MaterialPropertiesTable* GlassPropertiesTable = new G4MaterialPropertiesTable();
	GlassPropertiesTable->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexGlass, 2);
	GlassPropertiesTable->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLengthGlass, 2);
	BoronSilicateGlass->SetMaterialPropertiesTable(GlassPropertiesTable);

	//Air optical properties
	G4double RefractiveIndexAir[2] = {1.0002926, 1.0002926};
	G4double AbsorptionLengthAir[2] = {100. * m, 100. * m};

	G4MaterialPropertiesTable* AirPropertiesTable = new G4MaterialPropertiesTable();
	AirPropertiesTable->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexAir, 2);
	AirPropertiesTable->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLengthAir, 2);
	Air->SetMaterialPropertiesTable(AirPropertiesTable);

	/*	CHERENKOV WATER TANK	*/
	G4bool checkOverlaps = true;

	//World
	G4double world_sizeX = 7 * m;
	G4double world_sizeY = 7 * m;
	G4double world_sizeZ = 7 * m;

	G4Box* solidWorld = new G4Box("World_s", 0.5 * world_sizeX, 0.5 * world_sizeY, 0.5 * world_sizeZ);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World_l");
	G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0, checkOverlaps);

	//Water volume dimensions
	G4double WaterInnerRad = 0. * m, WaterOuterRad = 0.6 * m;
	G4double WaterHeight = 1.2 * m;

	//Polyethylene film dimensions
	G4double PEFilmThickness = 0.2 * mm;
	G4double PEFilmInnerRad = 0. * m, PEFilmOuterRad = 0.6 * m + 0.5 * PEFilmThickness;
	G4double PEFilmHeight = WaterHeight + PEFilmThickness;

	//Concrete tank dimensions
	G4double TankThickness = 100 * mm;
	G4double TankInnerRad = 0. * m, TankOuterRad = PEFilmOuterRad + TankThickness;
	G4double TankHeight = PEFilmHeight;

	//Soil cylinder dimensions
	G4double CyLThikness = 2 * m;
	G4double CylInnerRad = 0. * m, CylOuterRad = TankOuterRad + CyLThikness;
	G4double CylHeight = TankHeight;

	//Soil mount dimensions
	G4double MountUpInnerRad = 0. * m, MountDownInnerRad = 0. * m;
	G4double MountUpOuterRad = TankOuterRad;
	G4double MountDownOuterRad = CylOuterRad;
	G4double MountHeight = 2. * m;
	G4double StartPhi = 0. * deg, StopPhi = 360. * deg;

	//PMT glass dimensions 
	G4double UpperPartHAX = 133. * mm, UpperPartHAY = 133. * mm, UpperPartHAZ = 133. * mm;
	G4double UpperPartCut = 82. * mm, UpperPartDeposit = 82. * mm; 

	G4double MidConeUpOuterRad = 121 * mm, MidConeDownOuterRad = 104.8 * mm;
	G4double MidConeDownInnerRad = 0. * mm, MidConeUpInnerRad = 0. * mm;
	G4double MidConeHeight = 14.5 * mm;

	G4double MidEllipsHAX = 133. * mm, MidEllipsHAY = 133. * mm, MidEllipsHAZ = 79. * mm;
	G4double MidEllipseUpperCut = 35.25 * mm, MidEllipseDownCut = 32.8 * mm, MidEllipseDeposit = 35.25 * mm;

	G4double UpperConeUpOuterRad = 43. * mm, UpperConeDownOuterRad = 119. * mm;
	G4double UpperConeUpInnerRad = 0. * mm, UpperConeDownInnerRad = 0. * mm;
	G4double UpperConeHeight = 58. * mm;

	G4double GlassThickness = 2. * mm, PhotoCathodeThickness = 0.2 * mm;
	G4double GlassThicknessCone = 2.638 * mm;

	//Coordinates
	G4double XCyl = 0. * m, YCyl = 0. * m, ZCyl = 0. * m;																			//Soil cylinder coordinates
	G4double XMount = XCyl, YMount = YCyl, ZMount = 0.5 * (CylHeight + MountHeight);												//Soil mount coordinates
	G4double XTank = XCyl, YTank = YCyl, ZTank = ZCyl;																				//Concrete tank coordinates
	G4double XPEFilm = XCyl, YPEFilm = YCyl, ZPEFilm = ZCyl;																		//PE film coordinates
	G4double XWater = XCyl, YWater = YCyl, ZWater = ZCyl;																			//Water volume coordinates

	//PMT parts coordinates
	G4double XUpperCone = XCyl, YUpperCone = YCyl, ZUpperCone = 0.5 * (WaterHeight - UpperConeHeight);
	G4double XMidEllipse = XCyl, YMidEllipse = YCyl, ZMidEllipse = ZUpperCone - 0.5 * UpperConeHeight - MidEllipseDeposit;
	G4double XMidCone = XCyl, YMidCone = YCyl, ZMidCone = ZMidEllipse - MidEllipseDownCut - 0.5 * MidConeHeight;
	G4double XUpperPartGlass = XCyl, YUpperPartGlass = YCyl, ZUpperPartGlass = 0.5 * WaterHeight - UpperConeHeight - (MidEllipseUpperCut + MidEllipseDownCut) - MidConeHeight + UpperPartCut;

	//Volumes
	G4Cons* solidSoilMount = {nullptr}, * solidMidCone = {nullptr}, * solidUpperCone = {nullptr}, * solidMidConeVac = {nullptr}, * solidUpperConeVac = {nullptr}, * solidDownPhotCath = {nullptr};

	G4Ellipsoid* solidUpperPartGlass = {nullptr}, * solidMidEllipse = {nullptr}, * solidUpperPartVac = {nullptr}, * solidMidEllipseVac = {nullptr}, * solidUpPhotCath = {nullptr};

	G4Tubs* solidSoilCylinder = {nullptr}, * solidTank = {nullptr}, * solidWaterVolume = {nullptr}, * solidPEFilm = {nullptr};

	G4LogicalVolume* logicSoilMount = {nullptr}, * logicSoilCylinder = {nullptr}, * logicTank = {nullptr}, * logicWaterVolume = {nullptr}, * logicPEFilm = {nullptr}, 
		* logicUpperPartGlass = {nullptr}, * logicMidCone = {nullptr}, * logicMidEllipse = {nullptr}, * logicUpperCone = {nullptr}, * logicMidConeVac = {nullptr},
		 * logicUpperConeVac = {nullptr}, * logicUpperPartVac = {nullptr}, * logicMidEllipseVac = {nullptr}, * logicDownPhotCath = {nullptr}, * logicUpPhotCath = {nullptr};

	G4VPhysicalVolume* physSoilMount = {nullptr}, * physSoilCylinder = {nullptr}, * physTank = {nullptr}, * physWaterVolume = {nullptr}, * physPEFilm = {nullptr}, 
		* physUpperPartGlass = {nullptr}, * physMidCone = {nullptr}, * physMidEllipse = {nullptr}, * physUpperCone = {nullptr}, * physMidConeVac = {nullptr},
		 * physUpperConeVac = {nullptr}, * physUpperPartVac = {nullptr}, * physMidEllipseVac = {nullptr}, * physDownPhotCath = {nullptr}, * physUpPhotCath = {nullptr};

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

	//Building PMT
	solidUpperPartGlass = new G4Ellipsoid("UpperPartGlass_s", UpperPartHAX, UpperPartHAY, UpperPartHAZ, -UpperPartHAZ, - UpperPartCut);
	logicUpperPartGlass = new G4LogicalVolume(solidUpperPartGlass, BoronSilicateGlass, "UpperPartGlass_l");
	physUpperPartGlass = new G4PVPlacement(0, G4ThreeVector(XUpperPartGlass, YUpperPartGlass, ZUpperPartGlass), logicUpperPartGlass, "UPPER_PART_GLASS", logicWaterVolume, false, 0, checkOverlaps);

	solidMidCone = new G4Cons("MidCone_s", MidConeDownInnerRad, MidConeDownOuterRad, MidConeUpInnerRad, MidConeUpOuterRad, 0.5 * MidConeHeight, StartPhi, StopPhi);
	logicMidCone = new G4LogicalVolume(solidMidCone, BoronSilicateGlass, "MidCone_l");
	physMidCone = new G4PVPlacement(0, G4ThreeVector(XMidCone, YMidCone, ZMidCone), logicMidCone, "CONICAL_GLASS_PART", logicWaterVolume, false, 0, checkOverlaps);

	solidMidEllipse = new G4Ellipsoid("MidEllipse_s", MidEllipsHAX, MidEllipsHAY, MidEllipsHAZ, -MidEllipseDownCut, MidEllipseUpperCut);
	logicMidEllipse = new G4LogicalVolume(solidMidEllipse, BoronSilicateGlass, "MidEllipse_l");
	physMidEllipse = new G4PVPlacement(0, G4ThreeVector(XMidEllipse, YMidEllipse, ZMidEllipse), logicMidEllipse, "MID_ELLIPTICAL_PART_GLASS", logicWaterVolume, false, 0, checkOverlaps);

	solidUpperCone = new G4Cons("UpperCone_s", UpperConeDownInnerRad, UpperConeDownOuterRad, UpperConeUpInnerRad, UpperConeUpOuterRad, 0.5 * UpperConeHeight, StartPhi, StopPhi);
	logicUpperCone = new G4LogicalVolume(solidUpperCone, BoronSilicateGlass, "UpperCone_l");
	physUpperCone = new G4PVPlacement(0, G4ThreeVector(XUpperCone, YUpperCone, ZUpperCone), logicUpperCone, "UPPER_CONICAL_GLASS_PART", logicWaterVolume, false, 0, checkOverlaps);

	//Photocathode and vacuum inside PMT
	solidUpPhotCath = new G4Ellipsoid("UpperPartPhotCath_s", UpperPartHAX - GlassThickness, UpperPartHAY - GlassThickness, UpperPartHAZ - GlassThickness, -UpperPartHAZ, - UpperPartCut);
	logicUpPhotCath = new G4LogicalVolume(solidUpPhotCath, AlMaterial, "UpperPartPhotCath_l");
	physUpPhotCath = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicUpPhotCath, "PHOTOCATHODE_1", logicUpperPartGlass, false, 0, checkOverlaps);

	solidUpperPartVac = new G4Ellipsoid("UpperPartVac_s", UpperPartHAX - GlassThickness - PhotoCathodeThickness, UpperPartHAY - GlassThickness - PhotoCathodeThickness, 
		UpperPartHAZ - GlassThickness - PhotoCathodeThickness, -UpperPartHAZ, - UpperPartCut);
		
	logicUpperPartVac = new G4LogicalVolume(solidUpperPartVac, Vacuum, "UpperPartVac_l");
	physUpperPartVac = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicUpperPartVac, "VACUUM_1", logicUpPhotCath, false, 0, checkOverlaps);

	solidDownPhotCath = new G4Cons("DownPartPhotCath_s", MidConeDownInnerRad, MidConeDownOuterRad - GlassThicknessCone, MidConeUpInnerRad, MidConeUpOuterRad - GlassThicknessCone,
		 0.5 * MidConeHeight, StartPhi, StopPhi);

	logicDownPhotCath = new G4LogicalVolume(solidDownPhotCath, AlMaterial, "DownPartPhotCath_l");
	physDownPhotCath = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicDownPhotCath, "PHOTOCATHODE_2", logicMidCone, false, 0, checkOverlaps);

	solidMidConeVac = new G4Cons("MidConeVac_s", MidConeDownInnerRad, MidConeDownOuterRad - GlassThicknessCone - PhotoCathodeThickness, MidConeUpInnerRad, 
		MidConeUpOuterRad - GlassThicknessCone - PhotoCathodeThickness, 0.5 * MidConeHeight, StartPhi, StopPhi);

	logicMidConeVac = new G4LogicalVolume(solidMidConeVac, Vacuum, "MidConeVac_s");
	physMidConeVac = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicMidConeVac, "VACUUM_2", logicDownPhotCath, false, 0, checkOverlaps);

	solidMidEllipseVac = new G4Ellipsoid("MidEllipseVac_s", MidEllipsHAX - GlassThickness, MidEllipsHAY - GlassThickness, MidEllipsHAZ - GlassThickness, -MidEllipseDownCut, MidEllipseUpperCut);	
	logicMidEllipseVac = new G4LogicalVolume(solidMidEllipseVac, Vacuum, "MidEllipseVac_l");
	physMidEllipseVac = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicMidEllipseVac, "VACUUM_3", logicMidEllipse, false, 0, checkOverlaps);

	solidUpperConeVac = new G4Cons("UpperConeVac_s", UpperConeDownInnerRad, UpperConeDownOuterRad - GlassThicknessCone, UpperConeUpInnerRad, UpperConeUpOuterRad - GlassThicknessCone,
		 0.5 * UpperConeHeight, StartPhi, StopPhi);

	logicUpperConeVac = new G4LogicalVolume(solidUpperConeVac, Vacuum, "UpperConeVac_l");
	physUpperConeVac = new G4PVPlacement(0, G4ThreeVector(0. * mm, 0. * mm, 0. * mm), logicUpperConeVac, "VACUUM_4", logicUpperCone, false, 0, checkOverlaps);



	/*	OPTICAL SURFACES	*/

	//Border water - tyvek: diffuse reflection
	G4double ReflectivityWaterTyvek[2] = {0.98, 0.98};
	
	G4OpticalSurface* OpSurfaceWaterTyvek = new G4OpticalSurface("SurfaceWaterTyvek");
	OpSurfaceWaterTyvek->SetType(dielectric_dielectric);
	OpSurfaceWaterTyvek->SetFinish(groundfrontpainted);
	OpSurfaceWaterTyvek->SetModel(unified);

	G4MaterialPropertiesTable* SurWaterTyvekPT = new G4MaterialPropertiesTable();
	SurWaterTyvekPT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityWaterTyvek, 2);
	OpSurfaceWaterTyvek->SetMaterialPropertiesTable(SurWaterTyvekPT);

	G4LogicalBorderSurface* TyvekSurface = {nullptr};
	TyvekSurface = new G4LogicalBorderSurface("WaterTyvekSurface", physWaterVolume, physPEFilm, OpSurfaceWaterTyvek);

	//Border polyethylene - concrete: absorption
	G4double ReflectivityPEConcrete[2] = {0., 0.};
	
	G4OpticalSurface* OpSurfacePEConcrete = new G4OpticalSurface("SurfacePEConcrete");
	OpSurfacePEConcrete->SetType(dielectric_dielectric);
	OpSurfacePEConcrete->SetFinish(groundfrontpainted);
	OpSurfacePEConcrete->SetModel(unified);

	G4MaterialPropertiesTable* SurPEConcretePT = new G4MaterialPropertiesTable();
	SurPEConcretePT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityPEConcrete, 2);
	OpSurfacePEConcrete->SetMaterialPropertiesTable(SurPEConcretePT);

	G4LogicalBorderSurface* PESurface = {nullptr};
	PESurface = new G4LogicalBorderSurface("PEConcreteSurface", physPEFilm, physTank, OpSurfacePEConcrete);

	//Border boronsilicate glass - photocathode: mirror reflection
	G4double ReflectivityPhotocathode[2] = {0.99, 0.99};
	
	G4OpticalSurface* OpSurfaceGlassPhotocathode = new G4OpticalSurface("SurfaceGlassPhotocathode");
	OpSurfaceGlassPhotocathode->SetType(dielectric_metal);
	OpSurfaceGlassPhotocathode->SetFinish(ground);
	OpSurfaceGlassPhotocathode->SetModel(glisur);

	G4MaterialPropertiesTable* SurGlassPhotocathode = new G4MaterialPropertiesTable();
	SurGlassPhotocathode->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityPhotocathode, 2);
	OpSurfaceGlassPhotocathode->SetMaterialPropertiesTable(SurGlassPhotocathode);

	G4LogicalBorderSurface* UpperPhotCathSurface = {nullptr}, * LowerUpperPhotCathSurface = {nullptr};
	UpperPhotCathSurface = new G4LogicalBorderSurface("UpperPhotocathodeSurface", physUpperPartGlass, physUpPhotCath, OpSurfaceGlassPhotocathode);
	LowerUpperPhotCathSurface= new G4LogicalBorderSurface("LowerPhotocathodeSurface", physMidCone, physDownPhotCath, OpSurfaceGlassPhotocathode);

	/*	VISUAL ATTRIBUTES	*/

	//Colours
	G4Colour grey(0.5, 0.5, 0.5);
	G4Colour black(0.0, 0.0, 0.0);
	G4Colour red(1.0, 0.0, 0.0);
	G4Colour green(0.0, 1.0, 0.0);
	G4Colour blue(0.0, 0.0, 1.0);
	G4Colour cyan(0.0, 1.0, 1.0);
	G4Colour magenta(1.0, 0.0, 1.0);
	G4Colour yellow(1.0, 1.0, 0.0);

	//Making world invisible
	auto UniverseVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	UniverseVisAtt->SetVisibility(true);
	UniverseVisAtt->SetForceWireframe(true);

	logicWorld->SetVisAttributes(UniverseVisAtt);
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Water visual attributes
	auto WaterVisAtt = new G4VisAttributes(blue);
	WaterVisAtt->SetVisibility(true);
	WaterVisAtt->SetForceSolid(true);

	logicWaterVolume->SetVisAttributes(WaterVisAtt);
	//logicTank->SetVisAttributes(G4VisAttributes::GetInvisible());

	//Boronsilicate glass visual attributes
	auto GlassVisAtt = new G4VisAttributes(green);
	GlassVisAtt->SetVisibility(true);
	GlassVisAtt->SetForceSolid(true);

	logicUpperPartGlass->SetVisAttributes(GlassVisAtt);
	logicMidCone->SetVisAttributes(GlassVisAtt);
	logicUpperCone->SetVisAttributes(GlassVisAtt);
	logicMidEllipse->SetVisAttributes(GlassVisAtt);

	//Photocathode visual attributes
	auto PhotCathVisAtt = new G4VisAttributes(red);
	PhotCathVisAtt->SetVisibility(true);
	PhotCathVisAtt->SetForceSolid(true);

	logicDownPhotCath->SetVisAttributes(PhotCathVisAtt);
	logicUpPhotCath->SetVisAttributes(PhotCathVisAtt);

	//Vacuum visual attributes
	auto VacuumVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	VacuumVisAtt->SetVisibility(true);
	VacuumVisAtt->SetForceSolid(true);

	logicUpperConeVac->SetVisAttributes(VacuumVisAtt);
	logicMidConeVac->SetVisAttributes(VacuumVisAtt);
	logicUpperConeVac->SetVisAttributes(VacuumVisAtt);
	logicMidEllipseVac->SetVisAttributes(VacuumVisAtt);

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......