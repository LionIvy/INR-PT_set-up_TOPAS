#ifndef contouredScatterer_hh
#define contouredScatterer_hh

#include "TsVGeometryComponent.hh"

class contouredScatterer : public TsVGeometryComponent
{
public:
    contouredScatterer(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~contouredScatterer();

	G4VPhysicalVolume* Construct();

  G4double minG4(G4double* data);
  G4double maxG4(G4double* data);
};

#endif
