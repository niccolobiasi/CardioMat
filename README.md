# CardioMat
A Matlab toolbox for cardiac electrophysiology simulations

## Introduction

CardioMat is a Matlab toolbox for cardiac electrophysiology simulation focused on patient-specific modeling. It is designed to allow easy and fast construction of electrophysiology cardiac digital twins from segmented anatomical images. At the moment, CardioMat focused on ventricular simulations. The complete CardioMat pipeline is composed by the following steps:
- Geometry import: cardiac geometry must be provided as a 3D binary image as produce by image segmentation procedures. vtk2vox.m and adas2vox.m functions allows to convert vtk 3d mesh and ADAS3D shells to 3D binary image.
- Labeling of cardiac surface (tagSurface.m):  allows the user to manually select the base of the ventricles adn automatically determines epicardial and endocardial surfaces.
- Computation of ventricular coordinates (computeUVC.m): computes ventricular coordinates (apicobasal, transmural, transventricular) based on surface labels and user-selection of the apex of the ventricles.
- Generation of fiber orientation (generateFibers.m): generates fiber orientation with a rule-based method starting from ventricular coordinates.
- Generation of Purkinje network (createPurkije.m): create a Purkinje network based on ventricular coordinates with a constrained optimization method.
-  Run simulation (runSimulation.m): runs monodomain simulations with the options provided by the users.

## Usage

Each function is provided with an help documentation explaining the usage of the function. 3 example scripts showing the use of CardioMat in different scenarios are available at .....  

## Requirements

The runSimulation.m function uses the gpuArray structure thus it requires the Parallel Computing Toolbox.
The vox2carp.m function requires iso2mesh (dowload from:  toolbox available in the Matlab path). vox2carp.m allows to generates  
