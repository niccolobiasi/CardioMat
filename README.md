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
- Run simulation (runSimulation.m): runs monodomain simulations completely on GPU with the options provided by the users. For a complete list of the user-selectable options see the runSimulation help documentation. This also include a list of the currently available ionic models.

A similar pipeline for atrial models is available. See the example_atria script for details (available at [https://doi.org/10.5281/zenodo.14699421](https://doi.org/10.5281/zenodo.14699420)).
Methods are described in details in the associated publication (https://doi.org/10.1016/j.compbiomed.2024.109529).
A whole-heart simulation pipeline is now also available. This include an atrioventricular model connecting atria and His-Purkinje system.
For executing whole heart simulations the runSimulation_wh.m function must be used. Atrial and ventricular segmentation must be provided separately.
Methodological details about whole-heart models are provided in the associated publication (https://doi.org/10.1007/s00366-025-02214-z)

## Usage

Each function is provided with an help documentation explaining the usage of the function. Some example scripts showing the use of CardioMat in different scenarios are available at [https://doi.org/10.5281/zenodo.14699421](https://doi.org/10.5281/zenodo.14699420).  
These include scripts showing ventricular, atrial, and whole-heart models construction and simulation. Additionally, a script executing the N-version benchmark problem is provided. 
Simulation results are saved into a binary file and can be visualized with plotFrame.m and VolShFrame.m functions. The former uses standard Matlab figures, the latter is based on volshow() Matlab command and uses a viewer3d.
For whole-heart simulations the functions plotFrame_wh.m and VolShFrame_wh.m must be used. 
Please, note that not all CardioMat functionalities were extensively tested. We are developing CardioMat based on our research needs and some new features may have not been rigorously verified yet.
Thus, if you experience any bug, error, or inconsistent results, please open an issue and/or send an email to niccolo.biasi@ing.unipi.it with a description of the issue.

## Requirements

CardioMat toolbox was tested on Matlab R2023b and R2024a, and it is supposed to work for this and future releases.
The runSimulation.m function uses the gpuArray structure thus it requires the **Parallel Computing Toolbox**.
createPurkinje.m, VolSh.m, VolShFrame.m, and all dependant functions requires the **Image Processing Toolbox**.
The vox2carp.m function requires iso2mesh toolbox available in the Matlab path (dowload from: https://github.com/fangq/iso2mesh). vox2carp.m allows to generates 3D meshes in openCARP (https://opencarp.org/) format, including fiber orientation and domain labels.


## License 

CardioMat is licensed under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. You can redistribute CardioMat and/or modify it under its terms. 

## Citation

To cite CardioMat, please cite the associated publication "Biasi, Niccolò, et al. "A Matlab Toolbox for cardiac electrophysiology simulations on patient-specific geometries." Computers in Biology and Medicine 185 (2025): 109529" available at https://doi.org/10.1016/j.compbiomed.2024.109529.
