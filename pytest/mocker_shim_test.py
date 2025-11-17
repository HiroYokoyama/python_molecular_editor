import builtins
from unittest import mock
class _MockerShim:
    MagicMock = mock.MagicMock
    mock_open = mock.mock_open
    class _PatchCallable:
        def __call__(self, target, *args, **kwargs):
            p = mock.patch(target, *args, **kwargs)
            return p.start()
        def object(self, target, attribute, *args, **kwargs):
            p = mock.patch.object(target, attribute, *args, **kwargs)
            return p.start()
    patch = _PatchCallable()

m = _MockerShim()

# Do NOT patch builtins.open on import; that breaks pytest's assertion rewrite
# and can cause other silent errors. Provide explicit APIs to enable/disable
# the shim if a test or developer intends to use it.
_active_patch = None

def enable():
    """Activate the shim to mock builtin open().

    Returns the mock object created by mock_open so callers can adjust its
    return values if needed. Callers should call `disable()` when finished.
    """
    global _active_patch
    if _active_patch is not None:
        # already active
        return _active_patch
    new_open = m.mock_open()
    _active_patch = mock.patch('builtins.open', new_open)
    _active_patch.start()
    return new_open


def disable():
    """Stop the shim if it was previously enabled."""
    global _active_patch
    if _active_patch is None:
        return
    try:
        _active_patch.stop()
    finally:
        _active_patch = None
