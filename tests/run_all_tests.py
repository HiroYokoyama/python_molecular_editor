import sys
import os
import re
import argparse
import subprocess
import time

# Robust path setup
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(BASE_DIR, "moleditpy", "src")
UNIT_DIR = os.path.join(BASE_DIR, "tests", "unit")
INTEGRATION_DIR = os.path.join(BASE_DIR, "tests", "integration")
GUI_DIR = os.path.join(BASE_DIR, "tests", "gui")

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
    import traceback

    traceback.print_exc()


# ── Cross-platform teardown crash detection ──────────────────────────
# Qt/VTK cleanup can trigger fatal signals AFTER all tests have passed.
# We detect these by matching the exit code AND confirming pytest reported
# zero failures.
#
# Windows exit codes (from GetExitCodeProcess):
#   0xC0000005 → access violation  (-1073741819 signed)
#   0x80010108 → RPC_E_DISCONNECTED (-2147417848 signed)
#   0xC000013A → STATUS_CONTROL_C_EXIT (-1073741558 signed)
#
# Unix exit codes (negative = killed by signal):
#   -11 → SIGSEGV    -6 → SIGABRT    -5 → SIGTRAP
_KNOWN_CRASH_CODES = frozenset({
    # Windows (signed and unsigned representations)
    -1073741819,   # access violation (0xC0000005)
    -2147417848,   # RPC_E_DISCONNECTED (0x80010108)
    -1073741558,   # STATUS_CONTROL_C_EXIT (0xC000013A)
    3221225477,    # unsigned 0xC0000005 on some Python versions
    2147549704,    # unsigned 0x80010108 on some Python versions
    # Unix signals (subprocess returns -signal_number)
    -11,           # SIGSEGV
    -6,            # SIGABRT
    -5,            # SIGTRAP
})

# Pattern to extract pytest summary line, e.g. "5 passed, 1 warning in 0.52s"
_PYTEST_PASSED_RE = re.compile(r"(\d+) passed")
_PYTEST_FAILED_RE = re.compile(r"(\d+) failed")


def _is_teardown_crash(returncode, combined_output):
    """Check if a non-zero exit is a benign teardown crash (any platform).

    Returns True if the exit code matches a known crash AND pytest
    reported at least 1 passed test with 0 failures.
    """
    if returncode in _KNOWN_CRASH_CODES:
        passed = _PYTEST_PASSED_RE.search(combined_output)
        failed = _PYTEST_FAILED_RE.search(combined_output)
        has_passed = passed and int(passed.group(1)) > 0
        has_failed = failed and int(failed.group(1)) > 0
        if has_passed and not has_failed:
            return True
    return False


def run_suite(name, path, env_vars=None, extra_args=None, enable_cov=True):
    """Run a test suite in a separate process for isolation."""
    print(
        f"\n>>> Running {name} Tests{' with Coverage' if enable_cov else ''}...",
        flush=True,
    )

    # Use shorter timeout in headless CI mode for faster feedback
    timeout = "30" if os.environ.get("MOLEDITPY_HEADLESS") == "1" else "60"

    # Enable coverage and append to produce a combined .coverage file
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        "-vv",
        f"--timeout={timeout}",
        path,
        "--tb=short",
    ]

    if enable_cov:
        cmd.extend(["--cov=moleditpy", "--cov-append", "--cov-report="])

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
        start_ts = time.time()
        print(f"  [START] {name}: {' '.join(cmd)}", flush=True)

        process = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
        )
        output_lines = []
        assert process.stdout is not None
        for line in process.stdout:
            output_lines.append(line)
            # Stream immediately so CI logs show live progress / last test before freeze
            print(line, end="", flush=True)

        returncode = process.wait()
        duration = time.time() - start_ts
        print(
            f"  [END] {name}: exit={returncode} elapsed={duration:.1f}s",
            flush=True,
        )

        combined_output = "".join(output_lines)
        if returncode != 0:
            if _is_teardown_crash(returncode, combined_output):
                print(
                    f"  [NOTE] {name}: Ignoring post-test teardown crash "
                    f"(exit code {returncode}). All tests passed.",
                    flush=True,
                )
                return 0
        return returncode
    except Exception as e:
        print(f"Error running {name} tests: {e}")
        return 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unified Test Suite Runner")
    parser.add_argument(
        "--headless",
        action="store_true",
        help="Run GUI tests in headless mode (no Windows)",
    )
    parser.add_argument("--unit", action="store_true", help="Run ONLY Unit tests")
    parser.add_argument(
        "--integration", action="store_true", help="Run ONLY Integration tests"
    )
    parser.add_argument("--gui", action="store_true", help="Run ONLY GUI tests")
    parser.add_argument(
        "--no-report", action="store_true", help="Skip reporting phase entirely"
    )
    parser.add_argument(
        "--no-cov",
        action="store_true",
        help="Disable coverage collection (useful for CI without pytest-cov)",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Generate reports without running tests",
    )
    parser.add_argument(
        "--catalog-only", action="store_true", help="Update ONLY the assertion catalog"
    )
    parser.add_argument(
        "--skip-catalog",
        action="store_true",
        help="Skip updating catalog during reporting",
    )

    # Capture all remaining arguments to pass to pytest
    args, extra_pytest_args = parser.parse_known_args()
    if extra_pytest_args and extra_pytest_args[0] == "--":
        extra_pytest_args = extra_pytest_args[1:]

    # Handle --catalog-only early
    if args.catalog_only:
        print(">>> Updating Assertion Catalog only...", flush=True)
        subprocess.run(
            [
                sys.executable,
                os.path.join(
                    BASE_DIR, "tests", "utils", "generate_assertion_catalog.py"
                ),
            ],
            cwd=BASE_DIR,
        )
        sys.exit(0)

    # Clean up old coverage data to ensure a fresh combined report
    if not args.no_cov:
        cov_file = os.path.join(BASE_DIR, ".coverage")
        if not args.report_only and os.path.exists(cov_file):
            try:
                os.remove(cov_file)
            except Exception as e:
                print(f"Warning: Could not remove old .coverage file: {e}")

    # Handle --report-only
    if args.report_only:
        print(">>> Generating Reports only (from existing .coverage)...", flush=True)
        subprocess.run(
            [
                sys.executable,
                os.path.join(BASE_DIR, "tests", "utils", "print_cov.py"),
                "--skip-run",
            ],
            cwd=BASE_DIR,
        )
        if not args.skip_catalog:
            subprocess.run(
                [
                    sys.executable,
                    os.path.join(
                        BASE_DIR, "tests", "utils", "generate_assertion_catalog.py"
                    ),
                ],
                cwd=BASE_DIR,
            )
        sys.exit(0)

    env_vars = {}
    if args.headless:
        print("Running in HEADLESS mode (MOLEDITPY_HEADLESS=1)", flush=True)
        env_vars["MOLEDITPY_HEADLESS"] = "1"
        env_vars["QT_QPA_PLATFORM"] = "offscreen"

    enable_cov = not args.no_cov
    print(
        f"Starting Unified Test Suite (Coverage: {'ENABLED' if enable_cov else 'DISABLED'})...",
        flush=True,
    )
    if extra_pytest_args:
        print(f"Extra pytest arguments: {extra_pytest_args}")

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
            ret_code = run_suite(
                name,
                path,
                env_vars=env_vars,
                extra_args=extra_pytest_args,
                enable_cov=enable_cov,
            )
            results[name] = "PASSED" if ret_code == 0 else "FAILED"
        except KeyboardInterrupt:
            print(f"\nInterrupted during {name} tests.")
            results[name] = "INTERRUPTED"
            break
        except Exception as e:
            print(f"Unexpected error running {name} tests: {e}")
            results[name] = "ERROR"

    # Final Summary
    print("\n" + "=" * 40)
    print("         TEST SUITE SUMMARY")
    print("=" * 40)
    all_passed = True
    for name, path in [
        ("UNIT", UNIT_DIR),
        ("INTEGRATION", INTEGRATION_DIR),
        ("GUI", GUI_DIR),
    ]:
        status = results.get(name, "SKIPPED")
        print(f" {name:<12}: {status}")
        if status != "PASSED" and status != "SKIPPED":
            all_passed = False
    print("=" * 40)

    if all_passed:
        print("\nALL requested tests passed successfully!")

        # Only report if not disabled AND coverage was enabled
        if not args.no_report and enable_cov:
            print(
                "\n>>> Generating Final Reports (Markdown Coverage + Assertion Catalog)...",
                flush=True,
            )
            # Call print_cov.py with --skip-run because we already collected coverage!
            subprocess.run(
                [
                    sys.executable,
                    os.path.join(BASE_DIR, "tests", "utils", "print_cov.py"),
                    "--skip-run",
                ],
                cwd=BASE_DIR,
            )
            if not args.skip_catalog:
                subprocess.run(
                    [
                        sys.executable,
                        os.path.join(
                            BASE_DIR, "tests", "utils", "generate_assertion_catalog.py"
                        ),
                    ],
                    cwd=BASE_DIR,
                )

        sys.exit(0)
    else:
        print("\nSome tests failed or were interrupted. Check the output above.")
        sys.exit(1)
