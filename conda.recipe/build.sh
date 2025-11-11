#!/bin/bash

set -e

# Install the Python package
${PYTHON} -m pip install . -vv --no-deps --no-build-isolation

# Set executable permissions for the main entry point scripts
chmod +x ${PREFIX}/bin/sqanti3 || true
chmod +x ${PREFIX}/bin/sqanti3-qc || true
chmod +x ${PREFIX}/bin/sqanti3-filter || true
chmod +x ${PREFIX}/bin/sqanti3-rescue || true
chmod +x ${PREFIX}/bin/sqanti3-reads || true

# Copy the main Python scripts to bin if they're not already there
if [ -f sqanti3.py ]; then
    cp sqanti3.py ${PREFIX}/bin/ || true
fi
if [ -f sqanti3_qc.py ]; then
    cp sqanti3_qc.py ${PREFIX}/bin/ || true
fi
if [ -f sqanti3_filter.py ]; then
    cp sqanti3_filter.py ${PREFIX}/bin/ || true
fi
if [ -f sqanti3_rescue.py ]; then
    cp sqanti3_rescue.py ${PREFIX}/bin/ || true
fi
if [ -f sqanti3_reads.py ]; then
    cp sqanti3_reads.py ${PREFIX}/bin/ || true
fi

# Copy additional data and configuration files
if [ -f sqanti3_config.yaml ]; then
    cp sqanti3_config.yaml ${PREFIX}/bin/ || true
fi

echo "SQANTI3 installation complete"
