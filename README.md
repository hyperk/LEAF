# :seedling: LEAF :seedling:

<ins>**L**</ins>ow <ins>**E**</ins>nergy <ins>**A**</ins>lgorithm <ins>**F**</ins>ramework


# DESCRIPTION

This algorithm is an alternative and simple LE fitter than can be used for HyperK and SuperK.

~~~~~~~~~~~~~~~~~~~~~~~~~
2020/02/09: LEAF was convert as a C++ class and can be included in your code.
~~~~~~~~~~~~~~~~~~~~~~~~~


# DOWNLOAD AND PRE-REQUISITES

To download this repository use :

```
git clone --single-branch --branch develop_lite https://github.com/hyperk/LEAF.git
```

Then make sure the following pre-requisites are installed.
1. ROOT v5r34 or superior (not tested for older versions, but might work).
2. WCSim version compatible with your ROOT version.
3. BONSAI installation (although LEAF can work without it).
4. HKAstroAnalysis class is private and can be downloaded by SK collaborators on sukap cluster (although LEAF can work without it).


# COMPATIBILITY

Versions listed below have been tested so far.

1. WCSim-hybrid version: for geometries "HyperK", "HyperK_mPMT", "HyperK_HybridmPMT", "HyperK_HybridmPMT10PC"
2. In general, with all WCSim-hybrid geometries using whether BoxandLine20inchHQE or PMT3inchR14374 PMTs.
3. With official HK WCSim: Ask G. Pronost.


# INSTALLATION

1. Go to libWCSIM/ repository and make softlinks to your WCSim installation.
```
$ cd libWCSIM/
$ ln -s /path/to/your/WCSim/include/ include
$ ln -s /path/to/your/WCSim/src/ src
$ ln -s /path/to/your/WCSim/lib/libWCSimRoot.so libWCSimRoot.so
$ ln -s /path/to/your/WCSim/lib/libWCSimRoot.so.1.12.xx libWCSimRoot.so.1.12.xx
```
2. Go to the cloned LEAF repository and source the script RunAtStart.sh after making sure you are sourcing the proper ROOT directory.
```
$ cd /path/to/the/cloned/repository/LEAF
$ source ./RunAtStart.sh
```
3. Run the script ./SetupDataModel.sh to define the DataModel (if you have hk-AstroAnalysis, you should setup the global variable)
4. Enter the leaf/ repository and `make clean; make`
5. Enter the example/ repository and `make clean; make`
6. One example of how to run the code is set in example: test_example.sh


# TUNING FILES

Tuning files are provided in the inputs/ repository (e.g. ./inputs/timePDF_Directionality_DRnew.root) and are made from 10 MeV electrons generated isotropically and uniformly in the tank.
Alternatively, you can generate your own PDFs by following the commands below.

1. Go to macros/ repository and compile AnalyzeWSHierarchy and ProducePDF with GNUMake.
```
$ cd ./macros
$ make AnalyzeWSHierarchy
$ make ProducePDF
```

2. Produce plots using AnalyzeWSHierarchy which reads out WCSim output. Optionally, you can index a range of events to read with -s (start index) and -e (end index).
```
$ ./AnalyzeWSHierarchy  -i wcrim.root  -o plots.root  -s 0  -e 1000  -v
```

3. Produce time PDF and angular PDF using ProducePDF: uses plots made using AnalyzeWSHierarchy and generate PDFs for LEAF.
```
$ ./ProducePDF -i plots.root -o PDF.root
```

PathesREADME.md to PDF files required as inputs for LEAF are hard-coded inside leaf/LeafSplines.cc, inside the LoadSplines() methods.
Make sure these pathes are consistant with the files you generate on your own.



