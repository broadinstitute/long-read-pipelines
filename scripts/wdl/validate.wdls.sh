#!/bin/bash

# we do not set eu here because the WDLs may break and we want to catch that then continue

if [[ "$#" -eq 1 ]]; then
  cd $1
fi

echo "VALIDATING WDL FILES IN $PWD"
(womtool --version || echo "I need womtool to run" ) && echo

for wdl in $(find . -type f -name "*.wdl"); do
  womtool validate "${wdl}" 2>/dev/null
  test $? -eq 0 || \
    (echo -e "==============================\n${wdl}"; womtool validate "${wdl}")
done
