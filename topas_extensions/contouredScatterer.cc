// Component for contouredScatterer
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

#include "contouredScatterer.hh"

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

contouredScatterer::contouredScatterer(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


contouredScatterer::~contouredScatterer()
{
}
// TOPAS_EXTENSIONS_DIR = /home/ivy/work/TOPAS_projects/INR-PT_set-up_TOPAS/topas_extensions

G4VPhysicalVolume* contouredScatterer::Construct()
{
	BeginConstruction();

  //================================================================================
  // Fetch geometry parametrs
  //================================================================================

  G4int NbOfScLayers  = fPm->GetIntegerParameter( GetFullParmName("NbOfScLayers") );
  G4double* ScThickness     = fPm->GetDoubleVector( GetFullParmName("ScLayersHL"), "Length");
  G4double* Radii           = fPm->GetDoubleVector( GetFullParmName("ScLayersRmax"), "Length");
	G4String ScMaterial       = fPm->GetStringParameter( GetFullParmName("ScMaterial") );

  G4double CompensatorBoxHLX  = fPm->GetDoubleParameter(GetFullParmName("DegHLX"), "Length");
  G4double CompensatorBoxHLY  = fPm->GetDoubleParameter(GetFullParmName("DegHLY"), "Length");
	G4double CompensatorBaseHLZ  = fPm->GetDoubleParameter(GetFullParmName("DegBaseHLZ"), "Length");
  //const G4double CompensatorBoxHLZ  = fPm->GetDoubleParameter(GetFullParmName("DegrHLZ"), "Length");
  G4int NbOfDegLayers  = fPm->GetIntegerParameter( GetFullParmName("NbOfDegLayers") );
  G4double* CompensatorThickness = fPm->GetDoubleVector( GetFullParmName("DegLayersHL"), "Length");
  G4double* CompensatorRadii     = fPm->GetDoubleVector( GetFullParmName("DegLayersRmax"), "Length");
  G4String DegMaterial = fPm->GetStringParameter( GetFullParmName("DegMaterial") );

// G4double dummyX  = fPm->GetDoubleParameter(GetFullParmName("HLX"), "Length");
// G4double dummyY  = fPm->GetDoubleParameter(GetFullParmName("HLY"), "Length");
// G4double dummyZ  = fPm->GetDoubleParameter(GetFullParmName("HLZ"), "Length");

  G4double CompensatorBoxHLZ  = CompensatorBaseHLZ;
  if (NbOfDegLayers > 0){
		for(int i=0; i < NbOfDegLayers ; ++i){
    CompensatorBoxHLZ += CompensatorThickness[i];
    }
  }
  G4double TotalThickness = CompensatorBoxHLZ;
  if (NbOfScLayers > 0){
    for(int i=0; i < NbOfScLayers ; ++i){
      TotalThickness += ScThickness[i];
    }
  }
//  G4double TotalThickness = 5 *cm;

  //================================================================================
	// Geometry setup
	//================================================================================
	// Whole Box
  G4double HLX, HLY, HLZ;
  HLZ = TotalThickness;
	//HLY = HLX = 5*cm;
  maxG4(CompensatorRadii, NbOfDegLayers) > CompensatorBoxHLY ?  CompensatorBoxHLY = maxG4(CompensatorRadii, NbOfDegLayers) + 1.*mm : false;
  maxG4(Radii, NbOfScLayers) > CompensatorBoxHLY ? HLY = maxG4(Radii, NbOfScLayers) : HLY = CompensatorBoxHLY;

  maxG4(CompensatorRadii, NbOfDegLayers) > CompensatorBoxHLX ?  CompensatorBoxHLX = maxG4(CompensatorRadii, NbOfDegLayers) + 1.*mm : false;
  maxG4(Radii, NbOfScLayers) > CompensatorBoxHLX ? HLX = maxG4(Radii, NbOfScLayers) : HLX = CompensatorBoxHLX;




	G4String envelopeMaterialName = fParentComponent->GetResolvedMaterialName();
	G4Box* svWholeBox = new G4Box("TotalBox", HLX, HLY, HLZ);
	fEnvelopeLog = CreateLogicalVolume("TotalBox", envelopeMaterialName, svWholeBox);
	fEnvelopeLog->SetVisAttributes(fPm->GetInvisible());
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

  //--------------------------------------------------------------------------------
  // Scatterer
  //--------------------------------------------------------------------------------
	G4VisAttributes* red = new G4VisAttributes( G4Colour(1., 0., 0.));
    red -> SetVisibility(true);
    red -> SetForceSolid(true);

  G4double Z0Position = -HLZ;
  G4double ZCPosition;
	int bodyNb=0;
  for(int i = 0 ; i < NbOfScLayers ; ++i){
    if( fPm->ParameterExists(GetFullParmName("PrintInformation"))
			&&
			fPm->GetBooleanParameter(GetFullParmName("PrintInformation"))) {
			G4cout<<"Scattering Layer: \"" << i <<"\" , Thickness: "<< 2*ScThickness[i]/mm << " (mm), Radius: "<< Radii[i]/mm
				  <<" mm" <<G4endl;
		}
   G4Tubs* sFoilN = new G4Tubs("ScatteringFoil", 0, Radii[i], ScThickness[i], 0, 360*deg);
   G4LogicalVolume* lFoilN = CreateLogicalVolume("ScatteringFoil", ScMaterial , sFoilN);
   ZCPosition = Z0Position + ScThickness[i];
   Z0Position+= 2*ScThickness[i];
   G4ThreeVector* pos    = new G4ThreeVector(0.0, 0.0, ZCPosition);
   CreatePhysicalVolume("ScatteringFoil", i, true, lFoilN, 0, pos, fEnvelopePhys);
	 lFoilN->SetVisAttributes(red);
	 bodyNb = i;
  }
	//--------------------------------------------------------------------------------
  // Compensator Box Base
  //--------------------------------------------------------------------------------
	G4VisAttributes* lblue = new G4VisAttributes( G4Colour(204, 255, 238));
    lblue -> SetVisibility(true);
    lblue -> SetForceSolid(true);
		G4VisAttributes* lblue2 = new G4VisAttributes( G4Colour(204, 255, 238));
	    lblue2 -> SetVisibility(true);
	    lblue2 -> SetForceWireframe(true);

 if( fPm->ParameterExists(    GetFullParmName("PrintInformation"))		&&
		 		fPm->GetBooleanParameter(GetFullParmName("PrintInformation"))) {
		 		G4cout<<"Degrader base thickness: "<< 2*CompensatorBaseHLZ /mm << " (mm)" <<G4endl;
	}
	G4Box* sDegBaseBox = new G4Box("sDegBox", HLX, HLY, CompensatorBaseHLZ);
  G4LogicalVolume* lDegBaseBox = CreateLogicalVolume("lDegBaseBox", DegMaterial , sDegBaseBox);
	ZCPosition = Z0Position + CompensatorBaseHLZ;
	Z0Position+= 2*CompensatorBaseHLZ;
	G4ThreeVector* pos    = new G4ThreeVector(0.0, 0.0, ZCPosition);
	bodyNb++;
  CreatePhysicalVolume("DegBaseBox", bodyNb , true, lDegBaseBox, 0, pos, fEnvelopePhys);
	lDegBaseBox->SetVisAttributes(lblue);

  //--------------------------------------------------------------------------------
  // Compensator Box
  //--------------------------------------------------------------------------------
  G4double AirTubeR = maxG4(CompensatorRadii, NbOfDegLayers);
	std::cout << "/* AirTubeR = */" << AirTubeR << '\n';
	CompensatorBoxHLZ-=CompensatorBaseHLZ;
	G4VSolid* sDegBox = new G4Box ("sDegBox", HLX, HLY, CompensatorBoxHLZ);
	G4VSolid* hole    = new G4Tubs("hole", 0, AirTubeR, CompensatorBoxHLZ, 0, 360*deg);

	G4SubtractionSolid* sHoledDegBox = new G4SubtractionSolid("sDegBox", sDegBox, hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
	G4LogicalVolume* lDegBox = CreateLogicalVolume("lDegBox", DegMaterial , sHoledDegBox);
	ZCPosition = Z0Position + CompensatorBoxHLZ;
	//Z0Position+= 2*CompensatorBoxHLZ;
	G4ThreeVector* position = new G4ThreeVector(0,0, ZCPosition);
	bodyNb++;
	CreatePhysicalVolume("DegBox", bodyNb , true, lDegBox, 0, position, fEnvelopePhys);
	lDegBox->SetVisAttributes(lblue);

  //--------------------------------------------------------------------------------
  // Compensator "Rings"
  //--------------------------------------------------------------------------------

		for(int i = 0 ; i < NbOfDegLayers ; ++i){
			bodyNb++;
			if( fPm->ParameterExists(    GetFullParmName("PrintInformation")) &&
				  fPm->GetBooleanParameter(GetFullParmName("PrintInformation"))) {
				G4cout<<"Degrading Layer: \"" << i <<"\" , Thickness: "<< 2*CompensatorThickness[i]/mm << " (mm), Radius: "<< CompensatorRadii[i]/mm
						<<" mm" <<G4endl;
			}
	 if (CompensatorRadii[i] != AirTubeR){
	 	 G4Tubs* sDegN = new G4Tubs("sDegradingLvl", CompensatorRadii[i], AirTubeR, CompensatorThickness[i], 0, 360*deg);
	 	 G4LogicalVolume* lDegN = CreateLogicalVolume("lDegradingLvl", DegMaterial, sDegN);
		 ZCPosition = Z0Position + CompensatorThickness[i];
		 Z0Position+= 2*CompensatorThickness[i];
		 G4ThreeVector* pos3    = new G4ThreeVector(0.0, 0.0, ZCPosition);
		 CreatePhysicalVolume("DegradingLvl", bodyNb, true, lDegN, 0, pos3, fEnvelopePhys);
		 lDegN->SetVisAttributes(lblue);
	  }
	 // else{
		//  ZCPosition = Z0Position + CompensatorThickness[i];
		//  Z0Position+= 2*CompensatorThickness[i];
	 // }
	}




	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}
G4double contouredScatterer::minG4(G4double* data, G4int arrSize){
  //int arrSize = sizeof(data)/sizeof(data[0]);
  G4double res = data[0];
  if (arrSize == 0 ) return res;
  for(int i=1; i<arrSize; ++i){
    if (data[i] < res){
      res = data[i];
    }
  }
  return res;
}
G4double contouredScatterer::maxG4(G4double* data, G4int arrSize){
  //int arrSize = sizeof(data)/sizeof(data[0]);
  G4double res = data[0];
  if (arrSize == 0 ) return res;
  for(G4int i=1; i<arrSize; ++i){
    if (data[i] > res){
      res = data[i];
    }
  }
  return res;
}
