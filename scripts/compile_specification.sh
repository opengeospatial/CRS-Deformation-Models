#!/bin/sh
scriptfile=`realpath $0`
scriptdir=`dirname $scriptfile`
cd $scriptdir/..
if [ "$0" = "/metanorma/scripts/`basename $0`" ]; then
    # gem install asciidoctor-mathematical
    metanorma compile --agree-to-terms -t ogc -x xml,pdf,html,doc products/specification/abstract-specification-deformation-model-functional-model.adoc
    exit
fi
docker run --rm -v "$(pwd)":/metanorma -v "$(pwd)/.fontist/fonts/":/config/fonts metanorma/mn /metanorma/scripts/`basename $0`
sudo chown -R `id -u`:`id -g` .
