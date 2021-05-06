#!/bin/sh
asciidoctor --backend docbook --out-file - functional-model-strawman-cc.adoc | pandoc --from docbook --to docx --output functional-model-strawman-cc.docx

