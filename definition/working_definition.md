# Deformation model definition (for coordinate operations):

A deformation model is an approximate model of the displacement of the earth’s surface within a defined spatial and temporal extent.  It predicts the location of a point fixed to the surface  in terms of a specific accessible coordinate system at an any time within the temporal extent as a displacement from a time independent reference coordinate for that point.  

__Notes:__

* This defines a deformation model in the context of coordinate operations, specifically time dependent coordinate transformations or point motion models.  Time dependent transformations which do not involve deformation are also excluded, as these are covered by rotations or 14 parameter Bursa Wolf transformations (which can also express a constant rate of scale change)

* It may be more accurate to call this a displacement model, as it calculates displacements.  Mathematically deformation can be derived from it as the derivative of displacement with respect to position. The term deformation model is well established however, and appropriate as the model is only required in situations where there is ground deformation.

* The deformation model is only defined at the surface of the earth.  It is a function of horizontal position and time, but not of vertical position.  Geophysical models may predict deformation within the crust, but within the context of coordinate operations this has very little practical value.
.  
* “Accessible coordinate system” means a coordinate system in which locations can be measured directly, typically an ITRF coordinate system.  

* The deformation model may define horizontal and/or vertical displacement.  The deformation model is defined as a function calculating a displacement, rather than say a function that calculates the measureable coordinate directly from the reference coordinate.  However it must be accompanied by a specification of the formulae used to apply the displacement to a reference coordinate

* The deformation model should define formal uncertainties of displacements calculated from the model.  It may also define the correlation between displacements calculated for different locations and times. 

* The deformation model may represent deformation from multiple sources such as plate tectonic motion, glacial isostatic adjustment, seismic and co-seismic movement as well as human activities such as water or steam extraction.

* Practically a coordinate transformation implemented by a deformation model must be spatially continuous and uniquely invertible.  This limits the accuracy with which deformation such as fault rupture or very local movements will be represented by the model.

* The definition implies a “reference” coordinate system realised by the coordinates of points when the displacements are zero.  Commonly this is defined as the location of a point in a accessible coordinate system at a specific epoch (eg ITRF2014 at epoch 2020.0).  This is not accessible, except at that epoch.

* The deformation model must be accompanied by sufficient metadata to assess its appropriateness for use in a particular application.  This could include attributes such the authority and history of the model, the spatial and temporal extents it defines and the magnitude and dimension of the predicted displacements.  


# Deformation Model Functional Model definition:  

The deformation model functional model is the set of data, metadata, and algorithms required to fully use a deformation model, independently of any specific encoding or format of the model.
