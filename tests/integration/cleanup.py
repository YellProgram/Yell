import os
import shutil
import glob

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    work_dirs = glob.glob(os.path.join(base_dir, "work_*"))
    
    if not work_dirs:
        print("No work directories found to clean up.")
        return

    for d in work_dirs:
        print(f"Removing {d}")
        try:
            shutil.rmtree(d)
        except OSError as e:
            print(f"  Error removing {d}: {e}")

    print("Cleanup complete.")

if __name__ == "__main__":
    main()
