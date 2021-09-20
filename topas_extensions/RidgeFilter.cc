// Component for RidgeFilter
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

#include "RidgeFilter.hh"
#include "G4BooleanSolid.hh"
#include "TsParameterManager.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include <math.h>

RidgeFilter::RidgeFilter(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


RidgeFilter::~RidgeFilter()
{
}
// TOPAS_EXTENSIONS_DIR = /home/ivy/work/TOPAS_projects/INR-PT_set-up_TOPAS/topas_extensions

G4VPhysicalVolume* RidgeFilter::Construct()
{
	BeginConstruction();

  //================================================================================
  // Fetch geometry parametrs
  //================================================================================

  G4int NumberOfRidges            = fPm->GetIntegerParameter( GetFullParmName("NumberOfRidges") );
  G4int NumberOfModules           = fPm->GetIntegerParameter( GetFullParmName("NumberOfModules") );
  G4double* HeightOfRidges        = fPm->GetDoubleVector    ( GetFullParmName("HeightOfRidges"), "Length");
  G4double LenghtOfRidges         = fPm->GetDoubleParameter ( GetFullParmName("LenghtOfRidges"), "Length");
  G4double* WidthOfRidges         = fPm->GetDoubleVector    ( GetFullParmName("WidthOfRidges"), "Length");
  G4double DistanceBetweenModules = fPm->GetDoubleParameter ( GetFullParmName("DistanceBetweenModules"), "Length");
  G4String RFMaterial             = fPm->GetStringParameter( GetFullParmName("Material") );

  //================================================================================
  // Geometry setup
  //================================================================================
	// Whole Box
    G4double BlockWidth = sumG4(WidthOfRidges, NumberOfRidges) + DistanceBetweenModules;
    G4double HLX, HLY, HLZ;
    HLZ = 0.5*maxG4(HeightOfRidges, NumberOfRidges);
    HLX = 0.5*BlockWidth*NumberOfModules;
    HLY = 0.5*LenghtOfRidges;

	G4String envelopeMaterialName = fParentComponent->GetResolvedMaterialName();
	G4Box* svWholeBox = new G4Box("TotalBox", HLX, HLY, HLZ);
	fEnvelopeLog = CreateLogicalVolume("TotalBox", envelopeMaterialName, svWholeBox);
	fEnvelopeLog->SetVisAttributes(fPm->GetInvisible());
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

  //--------------------------------------------------------------------------------
  // RidgeFilter
  //--------------------------------------------------------------------------------
	G4VisAttributes* lblue = new G4VisAttributes( G4Colour(204., 255., 238.));
    lblue -> SetVisibility(true);
    lblue -> SetForceSolid(true);
    int BNum=0;

    G4double BlockPosition = -HLX+DistanceBetweenModules;
    G4double ElemPositionX;
    G4double ElemPositionZ;
    G4ThreeVector* pos;
    G4cout<<"BlockWidth " << BlockWidth  <<" cm" <<G4endl;


    for (int i = 0; i < NumberOfModules; ++i){

        ElemPositionX = BlockPosition;
            for (int j=0; j < NumberOfRidges; ++j){
                ElemPositionX+=0.5*WidthOfRidges[j];
                ElemPositionZ=-HLZ+0.5*HeightOfRidges[j];
                G4cout<<"ElemPositionX " << ElemPositionX  <<" cm" <<"== ElemPositionZ " << ElemPositionZ  <<" cm" <<G4endl;
                G4Box* sRidge = new G4Box("RidgeF", WidthOfRidges[j]*0.5, LenghtOfRidges*0.5, HeightOfRidges[j]*0.5);
                G4LogicalVolume* lRidge = CreateLogicalVolume("RidgeF", RFMaterial, sRidge);
                lRidge -> SetVisAttributes(lblue);
                pos = new G4ThreeVector(ElemPositionX, 0, ElemPositionZ);
                CreatePhysicalVolume("RFStair", BNum, true, lRidge, 0, pos, fEnvelopePhys);
                ElemPositionX+=0.5*WidthOfRidges[j];
                BNum++;
            }
         BlockPosition += BlockWidth;
     }

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}
G4double RidgeFilter::minG4(G4double* data, G4int arrSize){
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
G4double RidgeFilter::maxG4(G4double* data, G4int arrSize){
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
G4double RidgeFilter::sumG4(G4double* data, G4int arrSize){
  //int arrSize = sizeof(data)/sizeof(data[0]);
  G4double res = 0;
  if (arrSize == 0 ) return res;
  for(G4int i=0; i<arrSize; ++i){
    res += data[i];
  }
  return res;
}
