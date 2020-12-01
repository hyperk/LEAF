#!/bin/bash
#Script is for GitHub management, to make DataModel available for non-HK members

if [ "$HK_ASTROANALYSIS_DIR" == "" ] || [ -z $HK_ASTROANALYSIS_DIR ]; then
	echo "HK_ASTROANALYSIS_DIR is not defined or doesn't exist ($HK_ASTROANALYSIS_DIR). We can't create DataModel-lite."
else 
	echo "HK_ASTROANALYSIS_DIR was found to be: $HK_ASTROANALYSIS_DIR . Create DataModel-lite"
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/Geometry.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/GeoTools.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootEventInfo.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootHitInfo.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootPMTInfo.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootRecoInfo.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootTriggerInfo.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/TimeDelta.* .
	rsync -av  $HK_ASTROANALYSIS_DIR/DataModelRoot/RootPosition.* .
	
	rm *.o
		
	echo "You should update manually Environments.h"
fi
	
