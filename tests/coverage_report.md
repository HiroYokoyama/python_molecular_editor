........................................................................ [ 82%]
...............                                                          [100%]
============================== warnings summary ===============================
<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type SwigPyPacked has no __module__ attribute

<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type SwigPyObject has no __module__ attribute

<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type swigvarlink has no __module__ attribute

tests/unit/test_plugin_manager.py::test_plugin_info_version_tuple
tests/unit/test_plugin_manager.py::test_plugin_info_version_tuple
  E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\plugin_manager.py:425: DeprecationWarning: ast.Str is deprecated and will be removed in Python 3.14; use ast.Constant instead
    elif hasattr(ast, 'Str') and isinstance(node.value, ast.Str): # Py3.7 and below

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
87 passed, 5 warnings in 0.99s
....                                                                     [100%]
============================== warnings summary ===============================
<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type SwigPyPacked has no __module__ attribute

<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type SwigPyObject has no __module__ attribute

<frozen importlib._bootstrap>:488
  <frozen importlib._bootstrap>:488: DeprecationWarning: builtin type swigvarlink has no __module__ attribute

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
4 passed, 3 warnings in 0.32s
...............................................+++++++++++++++++++++++++++++++++++ Timeout +++++++++++++++++++++++++++++++++++
~~~~~~~~~~~~~~~~~~~~~~~~~ Stack of MainThread (8508) ~~~~~~~~~~~~~~~~~~~~~~~~~~
  File "<frozen runpy>", line 198, in _run_module_as_main
  File "<frozen runpy>", line 88, in _run_code
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pytest\__main__.py", line 9, in <module>
    raise SystemExit(pytest.console_main())
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\config\__init__.py", line 221, in console_main
    code = main()
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\config\__init__.py", line 197, in main
    ret: ExitCode | int = config.hook.pytest_cmdline_main(config=config)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_hooks.py", line 512, in __call__
    return self._hookexec(self.name, self._hookimpls.copy(), kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_manager.py", line 120, in _hookexec
    return self._inner_hookexec(hook_name, methods, kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_callers.py", line 121, in _multicall
    res = hook_impl.function(*args)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\main.py", line 365, in pytest_cmdline_main
    return wrap_session(config, _main)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\main.py", line 318, in wrap_session
    session.exitstatus = doit(config, session) or 0
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\main.py", line 372, in _main
    config.hook.pytest_runtestloop(session=session)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_hooks.py", line 512, in __call__
    return self._hookexec(self.name, self._hookimpls.copy(), kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_manager.py", line 120, in _hookexec
    return self._inner_hookexec(hook_name, methods, kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_callers.py", line 121, in _multicall
    res = hook_impl.function(*args)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\main.py", line 396, in pytest_runtestloop
    item.config.hook.pytest_runtest_protocol(item=item, nextitem=nextitem)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_hooks.py", line 512, in __call__
    return self._hookexec(self.name, self._hookimpls.copy(), kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_manager.py", line 120, in _hookexec
    return self._inner_hookexec(hook_name, methods, kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_callers.py", line 121, in _multicall
    res = hook_impl.function(*args)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 118, in pytest_runtest_protocol
    runtestprotocol(item, nextitem=nextitem)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 131, in runtestprotocol
    rep = call_and_report(item, "setup", log)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 244, in call_and_report
    call = CallInfo.from_call(
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 353, in from_call
    result: TResult | None = func()
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 245, in <lambda>
    lambda: runtest_hook(item=item, **kwds),
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_hooks.py", line 512, in __call__
    return self._hookexec(self.name, self._hookimpls.copy(), kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_manager.py", line 120, in _hookexec
    return self._inner_hookexec(hook_name, methods, kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_callers.py", line 121, in _multicall
    res = hook_impl.function(*args)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 165, in pytest_runtest_setup
    item.session._setupstate.setup(item)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\runner.py", line 523, in setup
    col.setup()
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\python.py", line 1723, in setup
    self._request._fillfixtures()
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 707, in _fillfixtures
    item.funcargs[argname] = self.getfixturevalue(argname)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 539, in getfixturevalue
    fixturedef = self._get_active_fixturedef(argname)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 627, in _get_active_fixturedef
    fixturedef.execute(request=subrequest)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 1110, in execute
    result: FixtureValue = ihook.pytest_fixture_setup(
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_hooks.py", line 512, in __call__
    return self._hookexec(self.name, self._hookimpls.copy(), kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_manager.py", line 120, in _hookexec
    return self._inner_hookexec(hook_name, methods, kwargs, firstresult)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\pluggy\_callers.py", line 121, in _multicall
    res = hook_impl.function(*args)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 1202, in pytest_fixture_setup
    result = call_fixture_func(fixturefunc, request, kwargs)
  File "C:\Users\hiro2\AppData\Local\Programs\Python\Python313\Lib\site-packages\_pytest\fixtures.py", line 908, in call_fixture_func
    fixture_result = next(generator)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\tests\gui\conftest.py", line 1049, in window
    main_window.close()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 242, in closeEvent
    widget.close()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 242, in closeEvent
    widget.close()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 242, in closeEvent
    widget.close()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 266, in closeEvent
    event.accept()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 266, in closeEvent
    event.accept()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 266, in closeEvent
    event.accept()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 266, in closeEvent
    event.accept()
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 597, in closeEvent
    return self.main_window_ui_manager.closeEvent(event)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window.py", line 134, in <lambda>
    return lambda *a, **k: attr(self._host, *a, **k)
  File "E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\moleditpy\src\moleditpy\modules\main_window_ui_manager.py", line 266, in closeEvent
    event.accept()
+++++++++++++++++++++++++++++++++++ Timeout +++++++++++++++++++++++++++++++++++
Wrote JSON report to E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\tests\cov_combined.json
Wrote HTML report to E:\Research\Calculation\moleditpy\20260212_01-moleditpy-test\python_molecular_editor\tests\coverage_html\index.html
============================================================
STEP 1: Running Unit Tests with coverage...
============================================================

============================================================
STEP 2: Running Integration Tests with coverage...
============================================================

============================================================
STEP 3: Running GUI Tests with coverage (headless)...
============================================================

============================================================
STEP 4: Generating combined coverage report...
============================================================
# Combined Coverage Report

**Suites**: Unit + Integration + GUI (Headless)

- **Total Statements**: 4402
- **Total Covered**: 1543
- **Overall Coverage**: **35.05%**

## File Breakdown

| File | Stmts | Miss | Cover |
| :--- | :--- | :--- | :--- |
| moleditpy\src\moleditpy\__init__.py                     |      1 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\__init__.py             |     18 |      6 |   66.7% |
| moleditpy\src\moleditpy\modules\analysis_window.py      |    114 |     25 |   78.1% |
| moleditpy\src\moleditpy\modules\atom_item.py            |    259 |    222 |   14.3% |
| moleditpy\src\moleditpy\modules\bond_item.py            |    341 |    308 |    9.7% |
| moleditpy\src\moleditpy\modules\calculation_worker.py   |    517 |    307 |   40.6% |
| moleditpy\src\moleditpy\modules\constants.py            |     27 |      0 |  100.0% |
| moleditpy\src\moleditpy\modules\main_window_app_state.py |    450 |    330 |   26.7% |
| moleditpy\src\moleditpy\modules\main_window_edit_actions.py |    968 |    821 |   15.2% |
| moleditpy\src\moleditpy\modules\main_window_molecular_parsers.py |    676 |    449 |   33.6% |
| moleditpy\src\moleditpy\modules\main_window_project_io.py |    256 |    190 |   25.8% |
| moleditpy\src\moleditpy\modules\main_window_string_importers.py |    173 |     35 |   79.8% |
| moleditpy\src\moleditpy\modules\mirror_dialog.py        |     65 |      8 |   87.7% |
| moleditpy\src\moleditpy\modules\molecular_data.py       |    203 |     26 |   87.2% |
| moleditpy\src\moleditpy\modules\plugin_interface.py     |     58 |     32 |   44.8% |
| moleditpy\src\moleditpy\modules\plugin_manager.py       |    276 |    100 |   63.8% |
| **TOTAL** | **4402** | **2859** | **35.1%** |

## Test Suite Status
