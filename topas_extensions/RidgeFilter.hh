#ifndef RidgeFilter_hh
#define RidgeFilter_hh

#include "TsVGeometryComponent.hh"

class RidgeFilter : public TsVGeometryComponent
{
public:
    RidgeFilter(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~RidgeFilter();

	G4VPhysicalVolume* Construct();

  G4double minG4(G4double* data, G4int arrSize);
  G4double maxG4(G4double* data, G4int arrSize);
  G4double sumG4(G4double* data, G4int arrSize);
};

#endif
