#!/bin/sh
#
# Script to compile the abstract specification document using the metanorma/mn docker image.  Runs docker 
# image providing the base directory as a virtual directory, then runs this script again inside the docker
# image to compile the specification document.
#
# Note that the docker image does not run (or I have not yet figured how to run) as a non root user, so 
# for the moment the ownership of output files is fixed after the script has completed.
#
scriptfile="$(realpath "$0")"
scriptdir="$(dirname "$scriptfile")"
cd "$scriptdir/.."
if [ "$0" = "/metanorma/scripts/$(basename "$0")" ]; then
    # gem install asciidoctor-mathematical
    metanorma compile --agree-to-terms -t ogc -x xml,pdf,html,doc products/specification/abstract-specification-functional-model-for-crustal-deformation.adoc
    exit
fi
docker run --rm -v "$(pwd)":/metanorma -v "$(pwd)/.fontist/fonts/":/config/fonts metanorma/metanorma "/metanorma/scripts/$(basename "$0")"
sudo chown -R "$(id -u):$(id -g)" .
