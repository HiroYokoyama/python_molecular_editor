## Summary

<!-- What does this PR do? List the key changes in 2–4 bullet points. -->

-
-

## Type of change

<!-- Check all that apply. -->

- [ ] Bug fix
- [ ] New feature
- [ ] Refactor / code cleanup
- [ ] Documentation
- [ ] Tests
- [ ] CI / build

## Related issues

<!-- Link any related issues, e.g. Closes #123 -->

## Test plan

<!-- How did you verify the change works? -->

- [ ] Ran `MOLEDITPY_HEADLESS=1 QT_QPA_PLATFORM=offscreen python tests/run_all_tests.py --no-cov --no-report --unit --integration` — all pass
- [ ] Manually tested in the GUI (describe what you tested below)
- [ ] Added new unit tests

**Manual test notes:**
<!-- Optional: describe what you clicked/drew/exported to verify -->

## Checklist

- [ ] `python -m flake8 moleditpy/src/ --select=F` returns 0 issues
- [ ] No bare `except:` clauses added (use `except Exception as e:` and log the error)
- [ ] No new unused imports or variables
- [ ] Plugin API changes are backwards-compatible (or documented in Summary)
