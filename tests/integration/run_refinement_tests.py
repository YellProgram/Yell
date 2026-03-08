import os
import shutil
import subprocess
import glob
import sys
import argparse
import re

# Configuration
DEFAULT_YELL_NEW = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../cmake-build-release/yell"))
if not os.path.exists(DEFAULT_YELL_NEW):
    DEFAULT_YELL_NEW = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../cmake-build-debug/yell"))

DEFAULT_THRESHOLD = 1e-4

def parse_parameters(file_path):
    """
    Parses Scale and RefinableVariables from a Yell input file or refined_parameters.txt.
    Returns a dictionary {name: value}.
    """
    params = {}
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Parse Scale - handle optional ESD in ()
    scale_match = re.search(r'Scale\s+([0-9.eE+-]+)(?:\([0-9]+\))?', content)
    if scale_match:
        params['Scale'] = float(scale_match.group(1))
    
    # Parse RefinableVariables block
    block_match = re.search(r'RefinableVariables\s*\[(.*?)\]', content, re.DOTALL)
    if block_match:
        block_content = block_match.group(1)
        # Find all name=value; possibly with ESD in () and possibly with comments after ;
        var_matches = re.findall(r'(\w+)\s*=\s*([0-9.eE+-]+)(?:\([0-9]+\))?\s*;', block_content)
        for name, value in var_matches:
            params[name] = float(value)
            
    return params

def run_refinement_test(target_file, starting_file, yell_bin, threshold):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    test_name = os.path.basename(target_file).replace("_target.txt", "").replace("_target copy.txt", "")
    work_dir = os.path.join(base_dir, f"work_refine_{test_name}")
    
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    
    print(f"Testing refinement: {test_name}")
    
    # 1. Generate "experimental" data from target
    with open(target_file, 'r') as f:
        target_content = f.read()
    # Force Refine false for data generation
    if "Refine true" in target_content:
        target_content = target_content.replace("Refine true", "Refine false")
    elif "Refine" not in target_content:
        target_content = "Refine false\n" + target_content
    
    with open(os.path.join(work_dir, "model.txt"), 'w') as f:
        f.write(target_content)
        
    print(f"  Generating experimental data from {os.path.basename(target_file)}...")
    try:
        subprocess.run([yell_bin], cwd=work_dir, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Data generation failed for {test_name}")
        print(f"  Stdout: {e.stdout}")
        print(f"  Stderr: {e.stderr}")
        return False
    
    if not os.path.exists(os.path.join(work_dir, "model.h5")):
        print(f"  ERROR: model.h5 not produced during data generation")
        return False
    
    shutil.move(os.path.join(work_dir, "model.h5"), os.path.join(work_dir, "experiment.h5"))
    
    # 2. Run refinement from starting point
    shutil.copy(starting_file, os.path.join(work_dir, "model.txt"))
    print(f"  Running refinement from {os.path.basename(starting_file)}...")
    try:
        subprocess.run([yell_bin], cwd=work_dir, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Refinement failed for {test_name}")
        print(f"  Stdout: {e.stdout}")
        print(f"  Stderr: {e.stderr}")
        return False
    
    if not os.path.exists(os.path.join(work_dir, "refined_parameters.txt")):
        print(f"  ERROR: refined_parameters.txt not produced")
        return False
    
    # 3. Compare parameters
    target_params = parse_parameters(target_file)
    refined_params = parse_parameters(os.path.join(work_dir, "refined_parameters.txt"))
    
    success = True
    print("  Comparing parameters:")
    all_keys = set(target_params.keys()) | set(refined_params.keys())
    for key in sorted(all_keys):
        t_val = target_params.get(key)
        r_val = refined_params.get(key)
        
        if t_val is None:
            print(f"    {key}: Extra parameter in refined results: {r_val}")
            success = False
        elif r_val is None:
            print(f"    {key}: Missing parameter in refined results (expected {t_val})")
            success = False
        else:
            diff = abs(t_val - r_val)
            status = "OK" if diff <= threshold else "FAIL"
            print(f"    {key}: target={t_val:g}, refined={r_val:g}, diff={diff:g} [{status}]")
            if diff > threshold:
                success = False
                
    if success:
        print(f"  SUCCESS")
    else:
        print(f"  FAILURE")
        
    return success

def main():
    parser = argparse.ArgumentParser(description="Run refinement integration tests for Yell.")
    parser.add_argument("--bin", default=DEFAULT_YELL_NEW, help=f"Path to yell binary (default: {DEFAULT_YELL_NEW})")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD, help=f"Parameter difference threshold (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("--filter", help="Filter tests by name (substring match)")
    
    args = parser.parse_args()

    if not os.path.exists(args.bin):
        print(f"ERROR: Yell binary not found at {args.bin}")
        sys.exit(1)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    models_dir = os.path.join(base_dir, "models_to_refine")
    
    target_files = glob.glob(os.path.join(models_dir, "*_target*.txt"))
    
    if args.filter:
        target_files = [f for f in target_files if args.filter in os.path.basename(f)]
    
    if not target_files:
        if args.filter:
            print(f"No target files matching filter '{args.filter}' found in {models_dir}")
        else:
            print(f"No target files (*_target.txt) found in {models_dir}")
        return

    failed = []
    for target_file in target_files:
        # Find matching starting file
        if "_target.txt" in target_file:
            starting_file = target_file.replace("_target.txt", "_starting.txt")
        elif "_target copy.txt" in target_file:
            starting_file = target_file.replace("_target copy.txt", "_starting.txt")
        else:
            # Fallback for any other variation
            starting_file = target_file.replace("target", "starting")
        
        if not os.path.exists(starting_file):
            print(f"WARNING: Starting file {starting_file} not found for {target_file}. Skipping.")
            continue
            
        if not run_refinement_test(target_file, starting_file, args.bin, args.threshold):
            failed.append(os.path.basename(target_file))
            
    if failed:
        print(f"\nRefinement tests FAILED for: {', '.join(failed)}")
        sys.exit(1)
    else:
        print("\nAll refinement integration tests PASSED!")

if __name__ == "__main__":
    main()
