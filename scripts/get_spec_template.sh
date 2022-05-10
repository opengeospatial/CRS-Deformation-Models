#!/bin/sh
#
# Script to extract the abstract specification template from the metanorma docker image.  This is run once only -
# just included here for the record as to how it was done.
#
docker run --rm -v "$(pwd)":/metanorma metanorma/mn metanorma new -d abstract-specification-topic -t ogc -l https://github.com/metanorma/mn-templates-ogc temp/specification
sudo chown -R "$(id -u):$(id -g)" .
