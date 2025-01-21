<p align="center">
  <img width="300" src="https://github.com/niccolobiasi/CardioMat/blob/main/logo_cardiomat_low.png">
</p>

# CardioMat
A Matlab toolbox for cardiac electrophysiology simulations

## Introduction

CardioMat is a Matlab toolbox for cardiac electrophysiology simulation focused on patient-specific modeling. It is designed to allow easy and fast construction of electrophysiology cardiac digital twins from segmented anatomical images. The complete CardioMat pipeline for ventricular models is composed by the following steps:

- Geometry import: cardiac geometry must be provided as a 3D binary image as produce by image segmentation procedures. vtk2vox.m and adas2vox.m functions allows to convert vtk 3d mesh and ADAS3D shells to 3D binary image.
- Labeling of cardiac surface (tagSurface.m):  allows the user to manually select the base of the ventricles adn automatically determines epicardial and endocardial surfaces.
- Computation of ventricular coordinates (computeUVC.m): computes ventricular coordinates (apicobasal, transmural, transventricular) based on surface labels and user-selection of the apex of the ventricles.
- Generation of fiber orientation (generateFibers.m): generates fiber orientation with a rule-based method starting from ventricular coordinates.
- Generation of Purkinje network (createPurkinje.m): create a Purkinje network based on ventricular coordinates with a constrained optimization method.
-  Run simulation (runSimulation.m): runs monodomain simulations with the options provided by the users.

A similar pipeline for atrial models is available. See the example_atria script for details (available at https://doi.org/10.5281/zenodo.14699421).
Methods are described in details in the associated publication (https://doi.org/10.1016/j.compbiomed.2024.109529).

## Usage

Each function is provided with an help documentation explaining the usage of the function. 3 example scripts showing the use of CardioMat in different scenarios are available at https://doi.org/10.5281/zenodo.14699421.  
Simulation results are saved into a binary file and can be visualized with plotFrame.m and VolShFrame.m functions. The former uses standard Matlab figures, the latter is based on volshow() Matlab command and uses a viewer3d. 

## Requirements

CardioMat toolbox was tested on Matlab R2023b and R2024a, and it is supposed to work for this and future releases.
The runSimulation.m function uses the gpuArray structure thus it requires the **Parallel Computing Toolbox**.
createPurkinje.m, VolSh.m, and VolShFrame.m requires the **Image Processing Toolbox**.
The vox2carp.m function requires iso2mesh toolbox available in the Matlab path (dowload from: https://github.com/fangq/iso2mesh). vox2carp.m allows to generates 3D mesh in openCARP (https://opencarp.org/) format, including fiber orientation and domain labels.


## License 

CardioMat is licensed under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. You can redistribute CardioMat and/or modify it under its terms. 

## Citation

To cite CardioMat, please cite the associated publication "Biasi, Niccolò, et al. "A Matlab Toolbox for cardiac electrophysiology simulations on patient-specific geometries." Computers in Biology and Medicine 185 (2025): 109529" available at https://doi.org/10.1016/j.compbiomed.2024.109529.
