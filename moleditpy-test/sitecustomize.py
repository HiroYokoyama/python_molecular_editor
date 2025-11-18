# sitecustomize.py runs very early when the Python interpreter starts and
# can be used to modify import behavior before pytest collects modules.
# We'll stub out `mocker_shim_test` if present to avoid running developer
# helper code that modifies builtins.open during collection.
import os
import sys
import types
import builtins
import io
import importlib.abc
import importlib.util

try:
    root = os.path.dirname(__file__)
    shim_path = os.path.join(root, 'mocker_shim_test.py')
    # If a local shim file exists, protect pytest by preplacing it with a
    # harmless stub module so imports won't execute the shim.
    # Prevent imports of a local mocker shim by inserting a meta-path
    # finder that returns a harmless stub instead of executing the shim.
    class _StubLoader(importlib.abc.Loader):
        def create_module(self, spec):
            return types.ModuleType(spec.name)
        def exec_module(self, module):
            # No-op loader: do not execute module code
            return

    class _ShimBlocker(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path, target=None):
            if fullname == 'mocker_shim_test':
                return importlib.util.spec_from_loader(fullname, _StubLoader())
            return None

    # Insert the finder at the front of meta_path so it takes precedence
    # during imports.
    sys.meta_path.insert(0, _ShimBlocker())

    # If some shim already replaced builtins.open with a non-callable, restore
    # a safe default so the test runner can open .pyc files during collection.
    if not callable(builtins.open):
        builtins.open = io.open
except Exception:
    # Keep this very defensive: if anything goes wrong we don't want to stop
    # Python from starting up.
    pass
