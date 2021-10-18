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
};

#endif
