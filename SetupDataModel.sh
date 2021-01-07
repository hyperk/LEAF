#!/bin/bash

if [ ! -z ./leaf/DataModel ]; then
	echo "Remove previous link"
	rm ./leaf/DataModel
fi

if [ "$HK_ASTROANALYSIS_DIR" == "" ] || [ -z $HK_ASTROANALYSIS_DIR ]; then
	echo "HK_ASTROANALYSIS_DIR is not defined or doesn't exist ($HK_ASTROANALYSIS_DIR). Will use lite DataModel."
	ln -s ${PWD}/leaf/DataModel-lite ./leaf/DataModel
else 
	echo "HK_ASTROANALYSIS_DIR was found to be: $HK_ASTROANALYSIS_DIR its DataModel will be used"
	ln -s $HK_ASTROANALYSIS_DIR/DataModelRoot ./leaf/DataModel
fi


