# CardioMat
A Matlab toolbox for cardiac electrophysiology simulations


CardioMat is a Matlab toolbox for cardiac electrophysiology simulation focused on patient-specific modeling. It is designed to allow easy and fast construction of electrophysiology cardiac digital twins from segmented anatomical images. At the moment, CardioMat focused on ventricular simulations. The complete CardioMat pipeline is composed by the following steps:
- Geometry import: cardiac geometry must be provided as a 3D binary image as produce by image segmentation procedures. vtk2vox.m and adas2vox.m functions allows to convert vtk 3d mesh and ADAS3D shells to 3D binary image.
- Labeling of cardiac surface (tagSurface.m):  allows the user to manually select the base of the ventricles adn automatically determines epicardial and endocardial surfaces.
- Computation of ventricular coordinates (computeUVC.m): 
