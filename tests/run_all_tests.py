import sys
import os
import argparse
import subprocess

# Robust path setup
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC_DIR = os.path.join(BASE_DIR, 'moleditpy', 'src')
UNIT_DIR = os.path.join(BASE_DIR, 'tests', 'unit')
GUI_DIR = os.path.join(BASE_DIR, 'tests', 'gui')

def run_suite(name, path, env_vars=None, extra_args=None):
    """Run a test suite in a separate process for isolation."""
    print(f"\n>>> Running {name} Tests...", flush=True)
    
    cmd = [sys.executable, "-m", "pytest", "-vv", "--timeout=60", path, "--tb=short"]
    if extra_args:
        cmd.extend(extra_args)
    
    env = os.environ.copy()
    if env_vars:
        env.update(env_vars)
    
    # Force src to be at the FRONT of PYTHONPATH to override site-packages
    pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = f"{SRC_DIR}{os.pathsep}{pythonpath}".strip(os.pathsep)
    # Also set it as the first entry in the current process sys.path just in case
    if SRC_DIR not in sys.path:
        sys.path.insert(0, SRC_DIR)

    try:
        result = subprocess.run(cmd, env=env, check=False)
        return result.returncode
    except Exception as e:
        print(f"Error running {name} tests: {e}")
        return 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified Test Suite Runner")
    parser.add_argument("--headless", action="store_true", help="Run GUI tests in headless mode (no Windows)")
    args = parser.parse_args()

    env_vars = {}
    if args.headless:
        print("Running in HEADLESS mode (MOLEDITPY_HEADLESS=1)", flush=True)
        env_vars["MOLEDITPY_HEADLESS"] = "1"
        env_vars["QT_QPA_PLATFORM"] = "offscreen"

        print("Starting Unified Test Suite (Unit + GUI)...", flush=True)
    
        unit_res = run_suite("UNIT", UNIT_DIR, env_vars=env_vars)
    
        gui_res = run_suite(
        "GUI", 
        GUI_DIR, 
        env_vars=env_vars
    )
    
    if unit_res == 0 and gui_res == 0:
        print("\nALL tests passed successfully!")
        sys.exit(0)
    else:
        print("\nSome tests failed. Check the output above.")
        sys.exit(1)
