#!/bin/bash
#Script is for GitHub management, to make DataModel available for non-HK members

if [ "$HK_ASTROANALYSIS_DIR" == "" ] || [ -z $HK_ASTROANALYSIS_DIR ]; then
	echo "HK_ASTROANALYSIS_DIR is not defined or doesn't exist ($HK_ASTROANALYSIS_DIR). We can't create DataModel-lite."
else 
	echo "HK_ASTROANALYSIS_DIR was found to be: $HK_ASTROANALYSIS_DIR . Create DataModel-lite"
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/Geometry.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/GeoTools.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/RootEventInfo.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/RootHitInfo.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/RootPMTInfo.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/RootRecoInfo.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/RootTriggerInfo.* .
	cp $HK_ASTROANALYSIS_DIR/DataModelRoot/TimeDelta.* .
	
	echo "You should update manually Environments.h"
fi
	
