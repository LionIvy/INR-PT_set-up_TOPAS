// Component for WallWithAHole
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

#include "wallwithahole.hh"

#include "G4BooleanSolid.hh"

#include "TsParameterManager.hh"

#include "G4Box.hh"

#include "G4Cons.hh"

#include "G4Tubs.hh"

#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"

#include "G4UIcommand.hh"

#include "G4SystemOfUnits.hh"
#include <math.h>

wallwithahole::wallwithahole(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
			 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
}


wallwithahole::~wallwithahole()
{
}


G4VPhysicalVolume* wallwithahole::Construct()
{
	BeginConstruction();

	//Call dimensions of Wall
	G4double HLX = fPm->GetDoubleParameter(GetFullParmName("HLX"), "Length");
	G4double HLY = fPm->GetDoubleParameter(GetFullParmName("HLY"), "Length");
	G4double HLZ = fPm->GetDoubleParameter(GetFullParmName("HLZ"), "Length");
	// G4double TransX = fPm->GetDoubleParameter(GetFullParmName("TransX"), "Length");
	// G4double TransY = fPm->GetDoubleParameter(GetFullParmName("TransY"), "Length");
	// G4double TransZ = fPm->GetDoubleParameter(GetFullParmName("TransZ"), "Length");
	// G4double RotX = fPm->GetDoubleParameter(GetFullParmName("RotX"), "Angle");
    // G4double RotY = fPm->GetDoubleParameter(GetFullParmName("RotY"), "Angle");
	// G4double RotZ = fPm->GetDoubleParameter(GetFullParmName("RotZ"), "Angle");


	//Create solid cylinder using dimensions
	G4VSolid* wall = new G4Box(fName, HLX, HLY, HLZ);

	//Call dimensions of cone
	G4double Rmin = fPm->GetDoubleParameter(GetFullParmName("Rmin"), "Length");
	G4double Rmax = fPm->GetDoubleParameter(GetFullParmName("Rmax"), "Length");
    //G4double HoleHalfLenth = fPm->GetDoubleParameter(GetFullParmName("HoleHalfLenth"), "Length");
    G4double HoleRotY = fPm->GetDoubleParameter(GetFullParmName("HoleRotY"), "Angle");
   // G4double totalHLZCone = fPm->SetDoubleParameter(GetFullParmName("HLCone"), "Length");
	//Create solid cone using dimensions
    G4VSolid* hole = new G4Tubs(fName, Rmin, Rmax, HLZ/sin(HoleRotY), 0, 360*deg);
	//(fName, RminCyl, RmaxCyl, totalHLZCyl, 0, 360*deg);

	//Define rotation and shift of cylinder relative to cone
	G4double shiftZ = 0;
	G4RotationMatrix* AlignMatrix = new G4RotationMatrix;
    AlignMatrix -> rotateY(- HoleRotY);
	const G4ThreeVector Vector(0,0,shiftZ);

	//Create solid union of cone and cylinder and rotate and shift cylinder relative to cone
	G4VSolid* fEnvelopeSolid = new G4SubtractionSolid(fName, wall, hole, AlignMatrix, Vector); //the first solid called in G4SubtractionSolid is the one that is placed relative to its mother ()TransX, RotX, ...). The second solid is placed relative to this first solid.

	//Create logical volume
	G4LogicalVolume* fEnvelopeLog = CreateLogicalVolume(fEnvelopeSolid);

	//Create physical volume
	G4VPhysicalVolume* fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	InstantiateChildren(fEnvelopePhys);

	return fEnvelopePhys;
}
