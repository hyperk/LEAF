# LEAF
Low Energy Algoritm Framework

This algorithm is an alternative and simple LE fitter than can be used for HyperK and SuperK.

~~~~~~~~~~~~~~~~~~~~~~~~~
2020/02/09: LEAF was convert as a C++ class and can be included in your code.

New compilation method:

	source RunAtStart.sh
	./SetupDataModel.sh
	cd leaf/
	make clean; make
	
In order to use the class in your code look at example/
~~~~~~~~~~~~~~~~~~~~~~~~~

# Pre-requisite to use the code:
1. BONSAI installation.
2. ROOT v5r34 or superior (not tested for older versions, but might work).
3. HKAstroAnalysis class is private and can be downloaded by SK collaborators on sukap cluster.

# Compatibility tested so far:
1. WCSim-hybrid version: for geometries "HyperK", "HyperK_mPMT", "HyperK_HybridmPMT", "HyperK_HybridmPMT10PC"
2. In general, with all WCSim-hybrid geometries using whether BoxandLine20inchHQE or PMT3inchR14374 PMTs.
3. With official HK WCSim: Ask G. Pronost.

# How to:
1. Source RunAtStart.sh after you updated your ROOT directory.
2. Use the script ./SetupDataModel.sh to define the DataModel (if you have hk-AstroAnalysis, you should setup the global variable)
2. Enter the leaf/ repository and make clean;make
3. Enter the example repository and make clean;make
4. One example of how to run the code is set in example: test_example.sh
5. inputs PDF, input from WCSim can be downloaded on sukap cluster. Please untar them in the LEAF repository.
6. You can use shell scripts in shell/ in order to run the fitter or launch on batch.

# Useful scripts in ./macros and ./shell
You can compile with GNUMake like following in ./macros:
```
$ make ProducePDF
```
## Making tuning file (e.g. ./inputs/timePDF_Directionality_DRnew.root)
1. Produce plots by AnalyzeWSHierarchy: reads out WCSim output and makes plots.
2. Produce time PDF (and angular PDF) by ProducePDF: uses plots made by AnalyzeWSHierarchy and generate PDFs for LEAF.

```
$ AnalyzeWSHierarchy -f wcrim_hybrid.root -o plots.root
$ ProducePDF -f plots.root -o PDF.root
```

## Making generic plots
- LEAFOutputAnalysisHybrid_leafclass: read LEAF output to produce generic plots. If one uses the master branch for LEAF, please use LEAFOutputAnalysisHybrid_master

## Shell scripts for analysis of large files
You can refere shell scripts in ./shell in order to analyze many files.
They are not working with latest LEAF class and its examples. They are just example how to analyze.

- generateShellXX.c: this is a root macro which generates shell scripts to be submitted to sukap by LauncherXX.sh
  - generateShell.c
  - generateShell_analyzeWCSim.c
- LauncherXX.sh: this submits jobs to sukap
  - Launcher.sh
  - Launcher_analyzeWCSim.sh
- Merger_analyzeWCSim.sh: merge generated output by Launcher_analyzeWCSim.sh

