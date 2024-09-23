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
    	2.885 * eV, 2.954 * eV, 3.026 * eV, 3.102 * eV, 3.181 * eV, 3.265 * eV
	};

	std::vector<G4double> RefractiveIndexWater = {
		1.3435, 1.344,  1.3445, 1.345,  1.3455, 1.346,  1.3465, 1.347,
    	1.3475, 1.348,  1.3485, 1.3492, 1.35,   1.3505, 1.351,  1.3518,
    	1.3522, 1.3530, 1.3535, 1.354,  1.3545, 1.355,  1.3555, 1.356
	}; 

	std::vector<G4double> AbsorptionLengthWater = {
		3.448 * m,  4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m,
    	15.152 * m, 17.241 * m, 18.868 * m, 20.000 * m, 26.316 * m, 35.714 * m,
    	45.455 * m, 47.619 * m, 52.632 * m, 52.632 * m, 55.556 * m, 52.632 * m,
    	52.632 * m, 47.619 * m, 45.455 * m, 41.667 * m, 37.037 * m, 33.333 * m
	};


	G4MaterialPropertiesTable* WaterOptProperties = new G4MaterialPropertiesTable();
	WaterOptProperties->AddProperty("RINDEX", PhotonEnergyWater, RefractiveIndexWater, PhotonEnergyWater.size());
	WaterOptProperties->AddProperty("ABSLENGTH", PhotonEnergyWater, AbsorptionLengthWater, PhotonEnergyWater.size());

	//BoronSilicate Glass optical properties
	G4double PhotonEnergy[2] = {1.9 * eV, 3.3 * eV};
	G4double RefractiveIndexGlass[2] = {1.5, 1.5};
	G4double AbsorptionLengthGlass[2] = {5. * m, 5. * m};

	G4MaterialPropertiesTable* GlassPropertiesTable = new G4MaterialPropertiesTable();
	GlassPropertiesTable->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexGlass, 2);
	GlassPropertiesTable->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLengthGlass, 2);

	//Air optical properties
	G4double RefractiveIndexAir[2] = {1.0002926, 1.0002926};
	G4double AbsorptionLengthAir[2] = {100. * m, 100. * m};

	G4MaterialPropertiesTable* AirPropertiesTable = new G4MaterialPropertiesTable();
	AirPropertiesTable->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexAir, 2);
	AirPropertiesTable->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionLengthAir, 2);

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
	logicUpPhotCath = new G4LogicalVolume(solidUpPhotCath, AlMaterial, "UpperPartPhotCth_l");
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

	return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......