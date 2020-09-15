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