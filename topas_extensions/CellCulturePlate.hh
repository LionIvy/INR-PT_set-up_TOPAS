#ifndef CellCulturePlate_hh
#define CellCulturePlate_hh

#include "TsVGeometryComponent.hh"

class CellCulturePlate : public TsVGeometryComponent //тут была ошибка
{
public:
    CellCulturePlate(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
    ~CellCulturePlate();

	G4VPhysicalVolume* Construct();

private:
  G4bool  wellsAreSolid = false;
  G4int   numberOfPlates = 1;

  G4int* platesOrientationVector;
  G4int plateOrientation = 1;
};

#endif
