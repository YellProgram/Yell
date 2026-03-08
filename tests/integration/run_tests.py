import os
import shutil
import subprocess
import glob
import sys
import argparse

# Configuration
DEFAULT_YELL_OLD = "yell"  # Should be on path
# Try to find the new version in common build locations
DEFAULT_YELL_NEW = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../cmake-build-release/yell"))
if not os.path.exists(DEFAULT_YELL_NEW):
    DEFAULT_YELL_NEW = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../cmake-build-debug/yell"))

DEFAULT_THRESHOLD = 0.001

def run_test(model_file, yell_old, yell_new, threshold):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    model_name = os.path.splitext(os.path.basename(model_file))[0]
    work_dir = os.path.join(base_dir, f"work_{model_name}")
    
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    
    shutil.copy(model_file, os.path.join(work_dir, "model.txt"))
    
    print(f"Testing model: {model_name}")
    
    # Run old version
    print(f"  Running old version ({yell_old})...")
    try:
        # We use check=True to catch failures, and capture_output to keep stdout clean
        subprocess.run([yell_old], cwd=work_dir, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Old version failed for {model_name}")
        print(f"  Return code: {e.returncode}")
        print(f"  Stderr: {e.stderr.decode()}")
        return False
    except FileNotFoundError:
        print(f"  ERROR: Old version '{yell_old}' not found on path.")
        return False
    
    if not os.path.exists(os.path.join(work_dir, "model.h5")):
        print(f"  ERROR: model.h5 not produced by old version")
        return False
    
    # Prepare experiment.h5 for the new version
    shutil.copy(os.path.join(work_dir, "model.h5"), os.path.join(work_dir, "experiment.h5"))
    
    # Run new version
    print(f"  Running new version ({yell_new})...")
    try:
        result = subprocess.run([yell_new], cwd=work_dir, check=True, capture_output=True, text=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: New version failed for {model_name}")
        print(f"  Return code: {e.returncode}")
        print(f"  Stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"  ERROR: New version '{yell_new}' not found.")
        return False
    
    # Parse Rw
    rw = None
    for line in output.splitlines():
        if " Rw=" in line:
            try:
                # Extract the value after Rw=
                parts = line.split(" Rw=")
                if len(parts) > 1:
                    rw_str = parts[1].split()[0] # Take first word after Rw=
                    rw = float(rw_str)
                break
            except (ValueError, IndexError):
                pass
    
    if rw is None:
        print(f"  ERROR: Could not find Rw in output")
        # Print a bit of output to help debugging
        print("  Output snippet:")
        print("\n".join(output.splitlines()[-10:]))
        return False
    
    print(f"  Rw = {rw}")
    if rw > threshold:
        print(f"  FAILURE: Rw {rw} > threshold {threshold}")
        return False
    
    print(f"  SUCCESS")
    return True

def main():
    parser = argparse.ArgumentParser(description="Run integration tests for Yell.")
    parser.add_argument("--old", default=DEFAULT_YELL_OLD, help=f"Path to old yell (default: {DEFAULT_YELL_OLD})")
    parser.add_argument("--new", default=DEFAULT_YELL_NEW, help=f"Path to new yell (default: {DEFAULT_YELL_NEW})")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"R-factor threshold (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("--models", help="Specific model file or directory to test")
    
    args = parser.parse_args()

    if not os.path.exists(args.new):
        print(f"ERROR: New yell binary not found at {args.new}")
        print("Please build the project first.")
        sys.exit(1)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    models_dir = os.path.join(base_dir, "models")
    
    if args.models:
        if os.path.isfile(args.models):
            models = [args.models]
        else:
            models = glob.glob(os.path.join(args.models, "*.txt"))
    else:
        models = glob.glob(os.path.join(models_dir, "*.txt"))

    if not models:
        print(f"No model files (.txt) found.")
        print(f"Please place your model files in {models_dir} or use --models")
        return

    failed = []
    for model in models:
        if not run_test(model, args.old, args.new, args.threshold):
            failed.append(os.path.basename(model))
            
    if failed:
        print(f"\nTests FAILED for: {', '.join(failed)}")
        sys.exit(1)
    else:
        print("\nAll integration tests PASSED!")

if __name__ == "__main__":
    main()
