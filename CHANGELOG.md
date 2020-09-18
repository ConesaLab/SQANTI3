# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.0] - 2020-09-21
### Fixed
- rebased on upstream master

## [1.4.9] - 2020-09-18
### Fixed
- Add missing `junctions` argument to `sqanti3_RulesFilter.main()

### Changed
- Reformatted CHANGELOG.md


## [1.4.8] - 2020-09-17
### Changed
- sqanti3_RulesFilter now has an explicit required `junctions` argument, replacing
  the implicit requirement that previously existed (where it assumed you were
  running the script in the directory of sqanti3_qc's output) to make running
  sqanti3_RulesFilter in pipelines clearer.
- sqanti3_qc now calls IsoAnnotLite_SQ1 directly instead of running it as a 
  script
- Incorporate commits [5e52b85](https://github.com/ConesaLab/SQANTI3/commit/5e52b85fc62474557b45618e054b07dfef580eb1)
  and [c10d715](https://github.com/ConesaLab/SQANTI3/commit/c10d7159288f4e7525eb7572f22e7d1812624741)
  from upstream


## [1.4.7] - 2020-09-15
### Fixed
- IsoAnnotLite_SQ1.py now works properly.

### Changed
- Update README.md to reflect using Miniconda to install or using the Docker
container


## [1.4.6] - 2020-09-14
### Changed
- Transitional update to clean up IsoAnnotLite_SQ1 using classes and context
  managers until I can get the Pandas-based version working


## [1.4.5] - 2020-09-14
### Fixed
- Change how the path to the "utilities" folder is determined
  (use `sqanti3.__path__[0]` instead of `os.path.dirname(__file__)`)

### Changed
- Revert `{NEWLINE}` and `{TAB}` to `\n` and `\t` in f-strings that was only
  was only changed in the first place because of my misunderstanding of their
  inclusion in f-strings


## [1.4.4] - 2020-09-13
### Fixed
- Fix docker issues
  - Remove the `-e` flag from the pip installation of SQANTI3
  - Update the versions required in the conda environment definition
  - Add using [mamba](https://github.com/mamba-org/mamba) to handle
conda package installation

### Added
- Add CHANGELOG.md

[1.4.9]: https://github.com/milescsmith/SQANTI3/compare/1.4.8...1.4.9
[1.4.8]: https://github.com/milescsmith/SQANTI3/compare/1.4.7...1.4.8
[1.4.7]: https://github.com/milescsmith/SQANTI3/compare/1.4.6...1.4.7
[1.4.6]: https://github.com/milescsmith/SQANTI3/compare/1.4.5...1.4.6
[1.4.5]: https://github.com/milescsmith/SQANTI3/compare/1.4.4...1.4.5
