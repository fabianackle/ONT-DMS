#!/usr/bin/env bash

nextflow run main.nf \
    -profile conda \
    -params-file params.json \
    -resume \
    -with-report \
    -with-timeline
