# Commit Backlog

This file tracks the history of commits made by the AI agent, providing a high-level summary of the changes and the reasoning behind them.

## [2026-05-21] - Fix Critical Bugs in SQANTI-reads Table and Plot Utility
- **Branch**: `sqanti_reads_plots`
- **Goal**: Fix critical execution bugs and code quality concerns in the SQANTI-reads tables and plots generation script.
- **Summary**:
    - Fixed a string-safe NaN comparison issue by replacing an invalid `np.isnan()` call with `pd.isna()`.
    - Removed duplicate `'junction_category'` column reference from `jxn_cols` schema.
    - Resolved a crashed column-mapping step caused by a leaked loop variable (`category` vs `cat`).
    - Replaced abrupt `sys.exit()` in library code with a standard `ValueError` exception and removed unreachable dead code check.
    - Updated comparison style for `None` to conform with PEP 8 (`is None` / `is not None`).
