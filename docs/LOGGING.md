# Logging in MoleditPy

This document describes the main application's logging system — how it is configured, where output goes, and the correct pattern for contributors writing internal app code.

---

## How Logging Is Set Up

Logging is configured **once** at startup in [`main.py`](../moleditpy/src/moleditpy/main.py) by `setup_logging()`, which is called before the Qt event loop starts.

```
main() → setup_logging() → QApplication → MainWindow
```

### What `setup_logging()` does

1. **Reads user settings** from `~/.moleditpy/settings.json` (before Qt is running):
   - `"log_to_file"` (bool) — whether to also write to a rotating log file.
   - `"log_level_debug"` (bool) — whether to use `DEBUG` level instead of `INFO`.

2. **Calls `logging.basicConfig(..., stream=sys.stdout)`** on the root logger with level `INFO` (or `DEBUG`) and a full format string:
   ```
   %(asctime)s [%(levelname)s] %(name)s (%(pathname)s:%(lineno)d): %(message)s
   ```

3. **Optionally adds a `RotatingFileHandler`** to the root logger:
   - Path: `~/.moleditpy/moleditpy.log`
   - Max size: 1 MB, 3 backup files kept.

4. **Installs `sys.excepthook`** so unhandled Python exceptions are logged as `ERROR` instead of printing to stderr silently.

5. **Routes Qt messages** (`qInstallMessageHandler`) to Python logging. Known noisy Qt warnings (e.g. `"Retrying to obtain clipboard"`) are silently downgraded to `DEBUG`.

---

## User-Facing Settings

| `settings.json` key | Type | Default | Effect |
|---|---|---|---|
| `log_to_file` | bool | `false` | Writes to `~/.moleditpy/moleditpy.log` in addition to stdout |
| `log_level_debug` | bool | `false` | Sets root logger level to `DEBUG` (very verbose) |

Settings are read **before Qt starts** so the full startup sequence is captured.

---

## The Correct Pattern for Internal Code

All modules inside `moleditpy/` must use a **named logger** — never the root logger directly.

```python
import logging

log = logging.getLogger(__name__)

# Then use it anywhere in the module:
log.debug("Detailed trace: %s", some_var)
log.info("File loaded: %s", path)
log.warning("Unexpected state; continuing anyway")
log.error("Failed to parse mol block: %s", e)
log.exception("Unhandled error in callback")  # includes full traceback
```

Using `__name__` means log records carry the fully-qualified module path (e.g. `moleditpy.ui.io_logic`), which makes filtering and debugging straightforward.

> [!CAUTION]
> **Never call `logging.basicConfig()` inside a module.** The root logger is already configured by `setup_logging()` at startup. A second `basicConfig()` call is silently ignored in Python 3, but it is misleading and may cause hard-to-debug behaviour if the call order ever changes.

> [!WARNING]
> **Never use `print()` for diagnostic output in production code.** `print()` bypasses the logging system entirely, does not appear in the log file, and cannot be suppressed by log level. Use `log.debug(...)` instead.

---

## Log Levels — When to Use Which

| Level | Call | Use when |
|---|---|---|
| `DEBUG` | `log.debug(...)` | Detailed trace info useful only when debugging (off by default) |
| `INFO` | `log.info(...)` | Normal significant events (file loaded, plugin initialised, etc.) |
| `WARNING` | `log.warning(...)` | Unexpected but recoverable situations |
| `ERROR` | `log.error(...)` | Failures that prevented an operation from completing |
| `CRITICAL` | `log.critical(...)` | App-level failures (rarely needed in module code) |

Use `log.exception(msg)` inside an `except` block — it logs at `ERROR` level and automatically appends the current traceback.

---

## Log File Location

```
~/.moleditpy/moleditpy.log          ← current
~/.moleditpy/moleditpy.log.1        ← previous (up to .3)
```

Only created when `"log_to_file": true` is set in `settings.json`. The file rotates automatically when it reaches 1 MB.

---

## Qt Message Integration

Qt's own logging (`qDebug`, `qWarning`, etc.) is piped through `_qt_message_handler` into Python logging at the equivalent level. All Qt messages appear with the prefix `Qt:` in the log output so they can be filtered separately.

Known noisy patterns (currently just `"Retrying to obtain clipboard"`) are downgraded to `DEBUG` to avoid flooding the log.

To add a new downgrade pattern, edit the `_DOWNGRADED_QT_PATTERNS` tuple in [`main.py`](../moleditpy/src/moleditpy/main.py):

```python
_DOWNGRADED_QT_PATTERNS = (
    "Retrying to obtain clipboard",
    "Your new noisy pattern here",
)
```
