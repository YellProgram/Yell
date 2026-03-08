# Yell Integration Tests

This directory contains scripts for regression testing Yell against an "old" version of the program.

## How it works

1. The script `run_tests.py` iterates through all `.txt` files in the `models/` folder.
2. For each model:
   - It creates a temporary directory `work_<model_name>`.
   - It copies the model file to `work_<model_name>/model.txt`.
   - It runs the **old** version of `yell` (found on the system PATH by default).
   - It takes the generated `model.h5` and copies it to `experiment.h5`.
   - It runs the **new** version of `yell` (found in `cmake-build-debug/yell` or `cmake-build-release/yell`).
   - It parses the output for ` Rw=...` (the R-factor).
   - It compares the R-factor against a threshold (default: 0.001).
   - If the R-factor is too large, it reports a failure.

## Usage

### 1. Prepare Models
Place any number of `model.txt` files into the `models/` directory, naming them `<something>.txt`.

### 2. Run Tests
From the project root:
```bash
python3 tests/integration/run_tests.py
```

Options:
- `--old <path>`: Path to the old version of `yell` (default: `yell`).
- `--new <path>`: Path to the new version of `yell` (default: auto-detected in build folders).
- `--threshold <value>`: R-factor threshold (default: 0.001).
- `--models <path>`: Path to a specific model file or a directory containing models.

### 3. Cleanup
To remove all `work_*` directories:
```bash
python3 tests/integration/cleanup.py
```

## Requirements
- Python 3.x
- `yell` (old version) must be on your shell PATH or specified via `--old`.
- The new version of Yell must be built in `cmake-build-debug` or `cmake-build-release`.
