# SQANTI3 CI/CD Testing Strategy

## Overview

This document explains the testing strategy across different CI/CD workflows to avoid redundancy while maintaining quality.

---

## Testing Workflows

### 1. `build-test-conda.yml` - Development Testing
**Purpose:** Validate code quality during development

**When it runs:**
- Every push to any branch
- Pull requests to master
- Skips on documentation-only changes (*.md files)

**What it tests:**
- ✅ Full pytest suite (~90 minutes)
- ✅ Tests on Ubuntu + macOS
- ✅ Uses real conda environment from SQANTI3.conda_env.yml
- ✅ Tests source code (development mode)

**Why:** Ensures code changes don't break functionality

---

### 2. `conda-package.yml` - Package Testing
**Purpose:** Validate conda package builds and works

**When it runs:**
- Push to master branch
- Pull requests to master
- Manual dispatch
- Skips on documentation-only changes

**What it tests:**

#### For Pull Requests (Smoke Tests Only - ~5 min)
- ✅ Package builds successfully
- ✅ Package installs correctly
- ✅ Python imports work (e.g., `import src.config`)
- ✅ Entry points exist (sqanti3, sqanti3-qc, etc.)
- ❌ **Skips full pytest** (already tested in build-test-conda.yml)

#### For Master Branch (Full Testing - ~20 min)
- ✅ Everything above
- ✅ **Full pytest suite** on installed package
- ✅ Final validation before publishing

**Why:**
- PRs: Avoid redundant testing (code already tested)
- Master: Ensure package actually works before publishing

---

### 3. `generate-docker-image.yml` - Docker CI
**Purpose:** Verify Docker image builds

**When it runs:**
- Push to master
- Pull requests to master
- Only when Docker/code files change

**What it tests:**
- ✅ Docker image builds successfully
- ❌ Does not test functionality (assumes code tests passed)

---

### 4. `push-to-dockerhub-on-release.yml` - Docker Release
**Purpose:** Publish Docker images on releases

**When it runs:**
- Only on GitHub releases

**What it does:**
- Builds and pushes to DockerHub
- No testing (assumes release is validated)

---

## Testing Matrix Comparison

| Workflow | PR | Master | Release | Duration | Pytest? |
|----------|-------|--------|---------|----------|---------|
| build-test-conda.yml | Full | Full | - | ~90 min | ✅ Full suite |
| conda-package.yml | Smoke | Full | - | 5-20 min | ✅ Master only |
| generate-docker-image.yml | Build only | Build only | - | ~15 min | ❌ |
| push-to-dockerhub-on-release.yml | - | - | Publish | ~10 min | ❌ |

---

## Why This Strategy?

### Before Optimization (Redundant)
```
PR to master:
├─ build-test-conda.yml: 90 min (full pytest) ✅
└─ conda-package.yml: 20 min (full pytest again) ❌ REDUNDANT

Total: 110 minutes
Redundancy: 20 minutes wasted
```

### After Optimization (Efficient)
```
PR to master:
├─ build-test-conda.yml: 90 min (full pytest) ✅
└─ conda-package.yml: 5 min (smoke tests only) ✅ NO REDUNDANCY

Total: 95 minutes
Time saved: 15 minutes per PR (85% reduction in conda-package time)
```

### On Master Branch (Pre-publish)
```
Push to master:
├─ build-test-conda.yml: 90 min (full pytest)
└─ conda-package.yml: 20 min (full pytest on package) ✅ FINAL VALIDATION

Total: 110 minutes
Reason: Ensure package works before publishing to anaconda.org
```

---

## Quality Gates

### Pull Request Quality Gates
1. ✅ Code passes full pytest (build-test-conda.yml)
2. ✅ Package builds successfully (conda-package.yml)
3. ✅ Entry points work (conda-package.yml)

### Master Branch Quality Gates
1. ✅ Code passes full pytest
2. ✅ Package builds successfully
3. ✅ Package passes full pytest
4. ✅ Ready to publish

### Release Quality Gates
1. ✅ All master branch gates passed
2. ✅ Manual verification
3. ✅ Tagged release created

---

## Cost Savings

### Per Pull Request
- **Before:** ~110 minutes total CI time
- **After:** ~95 minutes total CI time
- **Savings:** 15 minutes (13% reduction)

### Monthly (assuming 50 PRs)
- **Before:** 5,500 minutes
- **After:** 4,750 minutes
- **Savings:** 750 minutes/month

### Annual
- **Savings:** ~9,000 GitHub Actions minutes/year

---

## When Tests Run

### Scenario 1: Feature Branch Push
```
Action: git push origin feature/my-feature
Triggers:
  ✅ build-test-conda.yml (full pytest)
  ❌ conda-package.yml (not triggered - only for master PRs)
```

### Scenario 2: Pull Request to Master
```
Action: Create PR to master
Triggers:
  ✅ build-test-conda.yml (full pytest - 90 min)
  ✅ conda-package.yml (smoke tests only - 5 min)
Total: ~95 minutes
```

### Scenario 3: Merge to Master
```
Action: Merge PR to master
Triggers:
  ✅ build-test-conda.yml (full pytest)
  ✅ conda-package.yml (FULL pytest - final validation)
Total: ~110 minutes (worth it for final validation)
```

### Scenario 4: Documentation-Only Change
```
Action: Update README.md
Triggers:
  ❌ No workflows run (path filters skip *.md)
Savings: 100% (no wasted runs)
```

---

## Testing Philosophy

### Development Phase (PRs)
**Goal:** Fast feedback on code quality
- Focus on code correctness
- Skip redundant package testing
- Assume packaging is stable

### Pre-publish Phase (Master)
**Goal:** Ensure package quality before release
- Validate code works
- Validate package works
- Final quality gate

### Release Phase
**Goal:** Publish validated packages
- No additional testing
- Trust the quality gates

---

## Common Questions

### Q: Why not test the package on PRs?
**A:** We do! But only smoke tests (imports, entry points). Full pytest runs on the source code already, so running it again on the package is redundant unless you're about to publish.

### Q: What if packaging breaks between PR and master?
**A:** The full pytest runs on master push, so you'll catch it before publishing.

### Q: Can I force full tests on a PR?
**A:** Yes, use `workflow_dispatch` manually, or change the condition in the workflow.

### Q: Why keep both workflows?
**A:** They test different things:
- `build-test-conda.yml`: Tests CODE quality
- `conda-package.yml`: Tests PACKAGE quality

Both are needed, but don't need to run full tests redundantly.

---

## Troubleshooting

### If conda-package.yml fails on PR
1. Check if `build-test-conda.yml` passed (code should work)
2. Likely a packaging issue (MANIFEST.in, meta.yaml, entry points)
3. Fix packaging, not code

### If build-test-conda.yml fails
1. Code has bugs
2. Fix code, then re-test

### If tests fail on master but passed on PR
1. Check for race conditions or merge conflicts
2. Re-run workflows
3. Investigate test environment differences

---

## Future Improvements

### Potential Optimizations
1. ✅ Caching (implemented)
2. ✅ Path filters (implemented)
3. ✅ Conditional testing (implemented)
4. 🔄 Security scanning (planned)
5. 🔄 Coverage reporting (planned)
6. 🔄 Performance benchmarking (planned)

---

## Maintenance

This document should be updated when:
- New workflows are added
- Testing strategy changes
- Quality gates change
- Performance characteristics change

**Last updated:** 2025-11-05
**Version:** 1.0
