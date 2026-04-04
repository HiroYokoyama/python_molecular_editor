# GUI Testing Tips & Troubleshooting

This document contains best practices and solutions for common issues encountered during GUI testing with PyQt6 in MoleditPy, especially in headless environments.

## Infinite Recursion in CloseEvent (Headless Mode)

In headless environments (e.g., CI, WSL without X11), Qt's internal mechanism for suppressing duplicate close calls may fail. Calling `widget.close()` while processing a close event can trigger another close event, resulting in an infinite recursion.

### Symptoms
- `RecursionError: maximum recursion depth exceeded`
- Segmentation faults during test teardown.
- Tests hanging indefinitely during cleanup.

### Solution 1: Accept the Event (Recommended)
Instead of calling `close()` again inside a `closeEvent` or `eventFilter` catching a `Close` event, simply accept the event to let the Qt framework handle the window destruction naturally.

```python
def handle_close_event(self, event):
    # Perform any custom logic here (e.g., save data, show prompts)
    event.accept()  # Approve the close request
    return True     # Stop further event propagation
```

### Solution 2: Temporarily Remove the Event Filter
If your logic strictly requires calling `widget.close()` programmatically (e.g., during aggressive cleanup in `conftest.py`), you must remove the event filter immediately before calling it to prevent the filter from catching its own close event.

```python
def handle_close_event(self, event, widget):
    widget.removeEventFilter(self)  # Detaching filter to prevent recursion
    widget.close()                  # Now safe to call programmatically
    return True
```

### Direct Overrides
If you are overriding the `closeEvent(self, event)` method directly inside a `QWidget` subclass (rather than using an `eventFilter`), the exact same rule applies: you must call `event.accept()` instead of `self.close()`.

---
*Note: This tip was added based on issues discovered during headless testing on Windows/CI environments.*
