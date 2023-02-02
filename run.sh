#!/bin/bash

# Run cromwell with config file that routes tasks to disBatch.
java -Dconfig.file=local.conf -jar \
/mnt/home/mlek/mlek/bin/cromwell-56.jar run \
MultiSampleReadViz_consolidated.wdl \
-o workflow.options \
-i readviz_inputs_consolidated.json \
&> wdl_local.log


