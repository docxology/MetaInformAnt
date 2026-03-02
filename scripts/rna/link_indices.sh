#!/bin/bash
for d in output/amalgkit/shared/genome/*; do
  species=$(basename "$d")
  # Convert to amalgkit config name format (e.g. Apis_mellifera -> amellifera)
  # Actually, easier to let the pipeline handle it or just symlink directly based on the yaml dirs
done
