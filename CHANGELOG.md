# 2020-09-17: 1.4.8

## Changes
* sqanti3_RulesFilter now has an explicit required `junctions` argument, replacing
  the implicit requirement that previously existed (where it assumed you were
  running the script in the directory of sqanti3_qc's output) to make running
  sqanti3_RulesFilter in pipelines clearer.
* sqanti3_qc now calls IsoAnnotLite_SQ1 directly instead of running it as a 
  script
* Incorporate commits [5e52b85](https://github.com/ConesaLab/SQANTI3/commit/5e52b85fc62474557b45618e054b07dfef580eb1)
  and [c10d715](https://github.com/ConesaLab/SQANTI3/commit/c10d7159288f4e7525eb7572f22e7d1812624741)
  from upstream

# 2020-09-15: 1.4.7

## Fixes
* IsoAnnotLite_SQ1.py now works properly.

## Changes
* Update README.md to reflect using Miniconda to install or using the Docker
container


# 2020-09-14: 1.4.6

## Changes
* Transitional update to clean up IsoAnnotLite_SQ1 using classes and context
  managers until I can get the Pandas-based version working

# 2020-09-14: 1.4.5

## Fixes
* Change how the path to the "utilities" folder is determined
  (use `sqanti3.__path__[0]` instead of `os.path.dirname(__file__)`)

## Changes
* Revert `{NEWLINE}` and `{TAB}` to `\n` and `\t` in f-strings that was only
  was only changed in the first place because of my misunderstanding of their
  inclusion in f-strings

# 2020-09-13: 1.4.4

## Fixes
* Fix docker issues
  * Remove the `-e` flag from the pip installation of SQANTI3
  * Update the versions required in the conda environment definition
  * Add using [mamba](https://github.com/mamba-org/mamba) to handle
conda package installation

## Additions
* Add CHANGELOG.md