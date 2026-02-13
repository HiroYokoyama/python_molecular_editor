import sys
import os
import argparse
import subprocess

# Robust path setup
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC_DIR = os.path.join(BASE_DIR, 'moleditpy', 'src')
UNIT_DIR = os.path.join(BASE_DIR, 'tests', 'unit')
INTEGRATION_DIR = os.path.join(BASE_DIR, 'tests', 'integration')
GUI_DIR = os.path.join(BASE_DIR, 'tests', 'gui')

# Avoid colorama/COM issues on Windows by disabling color globally for pytest
os.environ["PYTEST_ADDOPTS"] = os.environ.get("PYTEST_ADDOPTS", "") + " --color=no"
os.environ["NO_COLOR"] = "1"

# Proactively neutralize colorama in the parent process if present
try:
    import colorama
    # convert=False prevents colorama from wrapping stdout/stderr with COM collectors on Windows
    colorama.init(convert=False)
except ImportError:
    pass
except Exception:
    pass

def run_suite(name, path, env_vars=None, extra_args=None):
    """Run a test suite in a separate process for isolation."""
    print(f"\n>>> Running {name} Tests...", flush=True)
    
    # Use shorter timeout in headless CI mode for faster feedback
    timeout = "30" if os.environ.get("MOLEDITPY_HEADLESS") == "1" else "60"
    cmd = [sys.executable, "-m", "pytest", "-vv", f"--timeout={timeout}", path, "--tb=short"]
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
    parser.add_argument("--unit", action="store_true", help="Run ONLY Unit tests")
    parser.add_argument("--integration", action="store_true", help="Run ONLY Integration tests")
    parser.add_argument("--gui", action="store_true", help="Run ONLY GUI tests")
    parser.add_argument("--no-report", action="store_true", help="Skip coverage and documentation reporting")
    args = parser.parse_args()

    env_vars = {}
    if args.headless:
        print("Running in HEADLESS mode (MOLEDITPY_HEADLESS=1)", flush=True)
        env_vars["MOLEDITPY_HEADLESS"] = "1"
        env_vars["QT_QPA_PLATFORM"] = "offscreen"

    print("Starting Unified Test Suite (Unit + Integration + GUI)...", flush=True)
    
    results = {}
    
    # Define suites to run
    suites = []
    run_all = not (args.unit or args.integration or args.gui)
    
    if args.unit or run_all:
        suites.append(("UNIT", UNIT_DIR))
    if args.integration or run_all:
        suites.append(("INTEGRATION", INTEGRATION_DIR))
    if args.gui or run_all:
        suites.append(("GUI", GUI_DIR))

    for name, path in suites:
        try:
            ret_code = run_suite(name, path, env_vars=env_vars)
            results[name] = "PASSED" if ret_code == 0 else "FAILED"
        except KeyboardInterrupt:
            print(f"\nInterrupted during {name} tests.")
            results[name] = "INTERRUPTED"
            break
        except Exception as e:
            print(f"Unexpected error running {name} tests: {e}")
            results[name] = "ERROR"

    # Final Summary
    print("\n" + "="*40)
    print("         TEST SUITE SUMMARY")
    print("="*40)
    all_passed = True
    for name, path in [("UNIT", UNIT_DIR), ("INTEGRATION", INTEGRATION_DIR), ("GUI", GUI_DIR)]:
        status = results.get(name, "SKIPPED")
        print(f" {name:<12}: {status}")
        if status != "PASSED" and status != "SKIPPED":
            all_passed = False
    print("="*40)
    
    if all_passed:
        print("\nALL requested tests passed successfully!")
        
        if not args.no_report:
            print("\n>>> Generating Final Reports (Coverage + Assertion Catalog)...", flush=True)
            subprocess.run([sys.executable, os.path.join(BASE_DIR, "tests", "utils", "print_cov.py")], cwd=BASE_DIR)
            
        sys.exit(0)
    else:
        print("\nSome tests failed or were interrupted. Check the output above.")
        sys.exit(1)
