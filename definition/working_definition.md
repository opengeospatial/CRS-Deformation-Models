# Deformation model definition (for coordinate operations):

A deformation model is an approximate model of the displacement of the earth’s surface within a defined spatial and temporal extent.  It predicts the location of a point fixed to the surface at an arbitrary time and in terms of an accessible coordinate system as an offset from a time invariant reference coordinate for that point.  

__Notes:__

* To be used confidently the deformation model must be accompanied by metadata that can be used to assess its fitness for use.  In particular this should define a measure of uncertainty associated with values calculated from the model.

* “Accessible coordinate system” means a coordinate system in which locations can be measured directly, typically an ITRF coordinate system.  

* The definition also implies a “reference” coordinate system within which coordinates of “fixed” points do not change over time.  Commonly this is defined as the location of a point in a accessible coordinate system at a specific epoch.  This is not accessible, except at that epoch.  The reference coordinate system may be used as if it were an accessible coordinate system, particularly for lower accuracy or local usage.

* The deformation model is only practially defined at the surface of the earth.  It is a function of horizontal position and time, but not of vertical position.  Geophysical models may predict deformation within the crust, but within the context of coordinate operations this has no practical value.

* The deformation model defines an offset between the accessible and reference coordinates. This may include both horizontal and vertical offsets.  Using an offset (rather than say a function that calculates the measureable coordinate directly from the reference coordinate) is a historical artefact as deformation models have traditionally come from geophysical models of ground displacements. 

* The deformation model may also define formal uncertainties of offsets calculated from the model.  It may also provide an estimate of correlation between offsets calculated for different locations or different times. 

* The deformation model may represent deformation from multiple sources such as plate tectonic motion, glacial isostatic adjustment, seismic and co-seismic movement as well as human activities such as water or steam extraction.

* Within the context of coordinate transformation a deformation model may be implemented as either a transformation between coordinate systems, or as a point motion model.

* Practically the coordinate transformation implemented by a transformation must be spatially continuous and uniquely invertible.  This  limits the accuracy with which deformation such as fault rupture or very local movements will be represented by the model.

* More accurately this could be called a displacement model as it predicts displacements (rather than deformation which may be calculated from it).  The term deformation model is well established however, and appropriate as it is only required in situations where there is ground deformation.

# Deformation Model Functional Model definition:  

The deformation model functional model is the set of data, metadata, and algorithms required to fully use a deformation model, independently of any specific encoding or format of the model.
