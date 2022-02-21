CRS-Deformation-Models
======================

This repository holds the artefacts and workings of the "Deformation model functional model" project team established by the CRS DWG at the Montreal virtual meeting on 15 June 2020. 

The next virtual meeting of this project team will be at 4pm EDT (20:00 UTC) on 14 March 2022 (https://www4.gotomeeting.com/join/616024989).  Meetings happen every four weeks. Details for joining the meeting will be posted here once they are available.

<!--
2020-09-07  https://www4.gotomeeting.com/join/494053021
2020-10-05 https://www4.gotomeeting.com/join/516068053
Note: EDT ends Nov 1
2020-11-02  https://www4.gotomeeting.com/join/882484381

2020-11-30  https://www4.gotomeeting.com/join/270784501
2020-12-28  https://www4.gotomeeting.com/join/437263613
2021-01-25  https://www4.gotomeeting.com/join/150445909
2021-02-22  https://www4.gotomeeting.com/join/577891581
2021-03-15  https://www4.gotomeeting.com/join/118008085
2021-04-12  https://www4.gotomeeting.com/join/694181949 20:00 UTC
2021-05-10	https://www4.gotomeeting.com/join/381926869
2021-06-07	https://www4.gotomeeting.com/join/792144133
2021-07-05	https://www4.gotomeeting.com/join/517573469
2021-08-02	https://www4.gotomeeting.com/join/403879845
2021-08-30	https://www4.gotomeeting.com/join/938424573

2021-09-27	https://www4.gotomeeting.com/join/449738333
2021-10-25      https://www4.gotomeeting.com/join/114555821
2021-11-22      https://www4.gotomeeting.com/join/864247981 21:00 UTC
2021-12-20      https://www4.gotomeeting.com/join/207460949

2022-01-17 https://www4.gotomeeting.com/join/544346533
2022-02-14 https://www4.gotomeeting.com/join/606539877
2022-03-14 https://www4.gotomeeting.com/join/616024989 20:00 UTC
2022-04-11 https://www4.gotomeeting.com/join/463903845

GGXF 2022-01-10 https://www4.gotomeeting.com/join/178306341
GGXF 2022-01-31 https://www4.gotomeeting.com/join/426755581
GGXF 2022-02-28 https://www4.gotomeeting.com/join/950826069
GGXF 2022-03-28 https://www4.gotomeeting.com/join/477678237
GGXF 2022-04-25 https://www4.gotomeeting.com/join/237877533

-->

This work is closely aligned with the work of the [GGXF (Geodetic Gridded data exchange format)](https://github.com/opengeospatial/CRS-Gridded-Geodetic-data-eXchange-Format) team.

## Active discussions

The project team is now working on compiling the abstract specification document that will be the main product of this teams work.

The document is open for comment and editing as a Google document at 
https://docs.google.com/document/d/1wlcB3zQjXqMmV-P6IHcLjUo-yC7PLeut_u9hazTbiEE/edit?usp=sharing.  This link can be used to submit comments on the document (select some text and then click the + icon to add a comment on the text).  

Other work done by this team is held in this repository, in particular the asciidoc source code for the abstract specification is in the products/specification folder.  The [github issues](https://github.com/opengeospatial/CRS-Deformation-Models/issues?q=is%3Aissue) also contain a body of discussion.

As an initial investigation this team undertook a survey of deformation models either planned or in use. Responses can be viewed at https://docs.google.com/spreadsheets/d/13IdqZDj8x8gVl7OTk7BkkI2HxfW61ng0y0WdVReD0us.  While the survey is formally closed, the questionnaire is still at https://docs.google.com/forms/d/11PCSVojPPD062P96veEjuqSKjnRMCrU8TyIruWCbWt0.  New information or updates are always welcome. 

## Project overview

The purpose of this team is to define the functional model(s), methods and workflows for computing coordinate offsets resulting from seismic and certain geodynamic and anthropogenic processes.  This work will be used to support developing
a standard for the deformation model functional model and/or format.

The motivation for this work is to support unambiguous and timely communication of deformation models from producers, typically national geodetic agencies, to the geospatial and positioning communities that require them.  
This is becoming increasingly critical as our dependence on global positioning systems and our requirements for accurate positioning increase.

The terms of reference for this project team proposes the following steps for this work. 

* Defining the term "deformation functional model" (DFM);
* Establishing the use case(s) of a DFM;
* Defining the user needs for a DFM;
* Evaluating existing deformation functional models and methods;
* Designing a DFM fit for its intended purpose as defined by the project team; 
* Ensure the DFM addresses model uncertainty and model validity;
* Develop a strategy for promulgating the model as a standard for the geodetic community.
* Develop an encoding structure for the functional model (when grid, refer to GGXF).

The work will be conducted under the auspices of the OGC CRS DWG in close collaboration with the CRS DWG project team developing the "Gridded Geodetic Exchange Format (GGXF)" and with IAG (International Association of Geodesy) WG 1.3.1 on "Time-dependent transformations between reference frames in deforming regions". 

## Collaboration

This is a public repository - everything is visible to anyone coming to this
website.  If you wish to be an active contributor with write access then you 
will need [join github](https://github.com/join) if you haven't already.  Once
you are signed on to github please raise a [new issue](https://github.com/opengeospatial/CRS-Deformation-Models/issues/new) with a request to be added.  If you would like to be included on the 
[Project Team](https://github.com/opengeospatial/CRS-Deformation-Models/wiki/Project-team) page then include a brief biography and a photo.  This is not 
required - be aware that everything on the page is visible to the public -  but it is good to have faces for names!  

Also once you have a github id you are encouraged to click the "Watch" button at the top of this page so that you will be notified of postings in the issues log.  Ideally we can capture most of the discussion in the [issue logs](https://github.com/opengeospatial/CRS-Deformation-Models/issues) where they will be recorded and easily searchable.  We will create an issue for each major discussion topic.


