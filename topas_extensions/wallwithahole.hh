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

#ifndef wallwithahole_hh
#define wallwithahole_hh

#include "TsVGeometryComponent.hh"

class wallwithahole : public TsVGeometryComponent
{    
public:
    wallwithahole(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~wallwithahole();
	
	G4VPhysicalVolume* Construct();
};

#endif
