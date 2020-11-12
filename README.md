# LEAF
Low Energy Algoritm Framework

This algorithm is an alternative and simple LE fitter than can be used for HyperK and SuperK.

~~~~~~~~~~~~~~~~~~~~~~~~~
2020/02/09: LEAF was convert as a C++ class and can be included in your code.

New compilation method:

	source RunAtStart.sh
	cd leaf/
	make clean; make
	
In order to use the class in your code look at example/

~~~~~~~~~~~~~~~~~~~~~~~~~

1. Source RunAtStart.sh after you updated your ROOT directory.

2. Enter the leaf/ repository and make clean;make

3. Enter the example repository and make clean;make

4. One example of how to run the code is set in example: test_example.sh

5. inputs PDF, input from WCSim can be downloaded on sukap cluster. Please untar them in the LEAF repository.

6. You can use shell scripts in shell/ in order to run the fitter or launch on batch.

Remarks:

1. HKAstroAnalysis class is private and can be downloaded by SK collaborators on sukap cluster.