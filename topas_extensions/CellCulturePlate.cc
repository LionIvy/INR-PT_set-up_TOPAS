// Component for CellCulturePlate
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "CellCulturePlate.hh"

#include "G4BooleanSolid.hh"

#include "TsParameterManager.hh"

#include "G4Box.hh"

#include "G4SubtractionSolid.hh"

#include "G4Tubs.hh"

#include "G4VisAttributes.hh"

#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"

#include "G4UIcommand.hh"

#include "G4SystemOfUnits.hh"
#include <math.h>

CellCulturePlate::CellCulturePlate(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
	if (fPm->ParameterExists(GetFullParmName("wellsAreSolid"))){
		wellsAreSolid = fPm -> GetBooleanParameter(GetFullParmName("wellsAreSolid"));
	}
	if (fPm->ParameterExists(GetFullParmName("numberOfPlates"))){
		numberOfPlates		= fPm -> GetIntegerParameter( GetFullParmName("numberOfPlates") );
	}else{
		numberOfPlates = 1;
	}
	if (fPm->ParameterExists(GetFullParmName("plateOrientation"))){
		plateOrientation		= fPm -> GetIntegerParameter( GetFullParmName("plateOrientation") );
	}else{
		plateOrientation = 1;
	}

	if (fPm->ParameterExists(GetFullParmName("platesOrientationVector"))){
		platesOrientationVector		= fPm -> GetIntegerVector( GetFullParmName("platesOrientationVector") );
	}else{
		platesOrientationVector = new G4int[numberOfPlates];
		for(G4int plateNum = 0; plateNum < numberOfPlates; ++plateNum){
			platesOrientationVector[plateNum] = plateOrientation;
		}
	}


}


CellCulturePlate::~CellCulturePlate()
{
}
G4VPhysicalVolume* CellCulturePlate::Construct()
{
	BeginConstruction();

  //================================================================================
  // Fetch geometry parametrs
  //================================================================================


	G4double	plateWidth				= fPm -> GetDoubleParameter ( GetFullParmName("plateWidth"), "Length");
  G4double	plateLength				= fPm -> GetDoubleParameter ( GetFullParmName("plateLength"), "Length");
	G4double	plateHeigth				= fPm -> GetDoubleParameter ( GetFullParmName("plateHeigth"), "Length");
	G4double	plateThickness		= fPm -> GetDoubleParameter ( GetFullParmName("plateThickness"), "Length");
  G4String	CCPMaterial				= fPm -> GetStringParameter( GetFullParmName("Material") );

  G4int numberOfColumns				= fPm -> GetIntegerParameter( GetFullParmName("numberOfColumns") );
  G4int numberOfRows					= fPm -> GetIntegerParameter( GetFullParmName("numberOfRows") );
  G4double wellsHeight				= fPm -> GetDoubleParameter ( GetFullParmName("wellsHeight"), "Length");
  G4double wellsInnerDiametr	= fPm -> GetDoubleParameter ( GetFullParmName("wellsInnerDiametr"), "Length");
	G4double wellsOuterDiametr	= fPm -> GetDoubleParameter ( GetFullParmName("wellsOuterDiametr"), "Length");
  G4double gapBetweenWells		= fPm -> GetDoubleParameter ( GetFullParmName("gapBetweenWells"), "Length");
  G4double wellsBottomThickness = fPm -> GetDoubleParameter ( GetFullParmName("wellsBottomThickness"), "Length");




	//================================================================================
  // Colors
  //================================================================================
	G4VisAttributes* yellowWire= new G4VisAttributes( G4Colour(255., 211., 0.));
    yellowWire -> SetVisibility(true);
    yellowWire -> SetForceWireframe(true);
	G4VisAttributes* yellowSolid= new G4VisAttributes( G4Colour(255., 211., 0.));
	  yellowSolid -> SetVisibility(true);
	  yellowSolid -> SetForceSolid(true);

	G4VisAttributes* blueFrame= new G4VisAttributes( G4Colour(0., 0., 255.));
		blueFrame -> SetVisibility(true);
		blueFrame -> SetForceWireframe(true);

	G4VisAttributes* whiteWire= new G4VisAttributes( G4Colour(255., 255., 255.));
	  whiteWire -> SetVisibility(true);
	  whiteWire -> SetForceWireframe(true);




  //================================================================================
  // Geometry setup
  //================================================================================
	// Whole Box
	//int numberOfPlates = 1;
	if (numberOfPlates	== 0){
		numberOfPlates	= 1;
	}


    G4double HLX, HLY, HLZ;
    HLZ = 0.5 * plateHeigth * numberOfPlates;
    HLX = 0.5 * plateWidth;//*127.6*mm;
    HLY = 0.5 * plateLength; //85.4*mm;

	G4String envelopeMaterialName = fParentComponent->GetResolvedMaterialName();
	G4Box* svWholeBox = new G4Box("TotalBox", HLX, HLY, HLZ);
	fEnvelopeLog = CreateLogicalVolume("TotalBox", envelopeMaterialName, svWholeBox);
	//fEnvelopeLog->SetVisAttributes(fPm->GetInvisible());
	fEnvelopeLog -> SetVisAttributes(blueFrame);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);


  //--------------------------------------------------------------------------------
  // the Plate
  //--------------------------------------------------------------------------------
  G4Box* plateBox;
	G4LogicalVolume* plateBoxLog;
	G4PVPlacement* plateBoxPhys;

	G4Box* emptyBox;
	G4LogicalVolume* emptyBoxLog;
	G4PVPlacement* emptyBoxPhys;

  int wellNum = 0;
	G4Tubs* sWell;
	G4LogicalVolume* lWell;

	G4Tubs* sWellBottom;
	G4LogicalVolume* lWellBottom;

	G4double wellsAreaX = 2* (0.5 * plateWidth - plateThickness);
	G4double wellsAreaY = 2* (0.5 * plateLength - plateThickness);
	G4double wellsAreaZ = 2* (0.5 * plateHeigth - plateThickness);
//|
// |  wellsAreaX |
// |+0-0-0-0-0-0+|
// |dX0-0-0-0-0-0dX|
G4double dX, dY; //Расстояние от границы рабочей зоны до первой лунки
dX = 0.5 * (wellsAreaX - numberOfColumns * wellsOuterDiametr - (numberOfColumns - 1) * gapBetweenWells);
dY = 0.5 * (wellsAreaY - numberOfRows * wellsOuterDiametr - (numberOfRows - 1) * gapBetweenWells);

if(dX < 0){
	G4cout<<"collumns size mismatch, distance between plate volume and fist well = " << dX/mm  <<" cm" <<G4endl;
	return 0;
}
if(dY < 0){
	G4cout<<"rows size mismatch, distance between plate volume and fist well = " << dY/mm  <<" cm" <<G4endl;
	return 0;
}
G4double wellPositionX, wellPositionY;
G4double wellPositionZ = 0.5 * (wellsAreaZ - wellsHeight);
G4ThreeVector* pos;
plateBox = new G4Box("plateBox", 0.5 * plateWidth, 0.5 * plateLength, 0.5 * plateHeigth);
plateBoxLog = CreateLogicalVolume("logicPlateBox", CCPMaterial, plateBox);
plateBoxLog -> SetVisAttributes(yellowWire);

emptyBox = new G4Box("emptyBox", 0.5 * wellsAreaX, 0.5 * wellsAreaY, 0.5 * wellsAreaZ);
emptyBoxLog = CreateLogicalVolume("logicEmptyBox", envelopeMaterialName, emptyBox);
emptyBoxLog -> SetVisAttributes(whiteWire);

sWell = new G4Tubs("Well", wellsInnerDiametr*0.5, wellsOuterDiametr*0.5, wellsHeight*0.5, 0., 360*deg);
lWell = CreateLogicalVolume("Well", CCPMaterial, sWell);
if(wellsAreSolid){
			lWell -> SetVisAttributes(yellowSolid);}
else{ lWell -> SetVisAttributes(yellowWire); }

sWellBottom = new G4Tubs("WellBottom", 0, wellsInnerDiametr*0.5, wellsBottomThickness*0.5, 0., 360*deg);
lWellBottom = CreateLogicalVolume("Well", CCPMaterial, sWellBottom);
lWellBottom -> SetVisAttributes(yellowSolid);


G4RotationMatrix* rMatrix			= new G4RotationMatrix();
G4RotationMatrix* rMatrix180	= new G4RotationMatrix();
rMatrix180 -> rotateX(180.0*deg);


for (G4int plateNum = 0; plateNum < numberOfPlates; ++plateNum){

	if (numberOfPlates > 1){
		if(platesOrientationVector[plateNum] >= 0){
			rMatrix = 0;
		}else{
			rMatrix = rMatrix180;
		}
	}else{
		if(plateOrientation >= 0){
			rMatrix = 0;
		}else{
			rMatrix = rMatrix180;
		}
	}


	plateBoxPhys = new G4PVPlacement(rMatrix,G4ThreeVector(0, 0, -HLZ + 0.5 * plateHeigth + plateNum * plateHeigth),
																		"physicalPlateBox",
																		plateBoxLog,
																		fEnvelopePhys,false,0);

	emptyBoxPhys = new G4PVPlacement(0,
																					G4ThreeVector(0, 0, -0.5 * (plateHeigth-wellsAreaZ)),
																					"physicalEmptyBox",
																					emptyBoxLog,
																					plateBoxPhys,false,0);

	for (G4int col = 0; col < numberOfColumns; ++col){
	    wellPositionX =  -0.5 * wellsAreaX + dX + 0.5 * wellsOuterDiametr + (wellsOuterDiametr + gapBetweenWells) * col;
	        for (G4int row = 0; row < numberOfRows; ++row){
	            wellPositionY =  -0.5 * wellsAreaY + dY + 0.5 * wellsOuterDiametr + (wellsOuterDiametr + gapBetweenWells) * row;
							pos = new G4ThreeVector(wellPositionX, wellPositionY, wellPositionZ);
	            CreatePhysicalVolume("CCPWell", wellNum, true, lWell, 0, pos, emptyBoxPhys);
	            wellNum++;
							pos = new G4ThreeVector(wellPositionX, wellPositionY, wellPositionZ-wellsHeight*0.5+wellsBottomThickness*0.5);
							CreatePhysicalVolume("WellBottom", wellNum, true, lWellBottom, 0, pos, emptyBoxPhys);
	            wellNum++;
	        }
	 }
}

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}
