# 25. Python is 3.13 only

`requires-python` floored at 3.12 while the lock held one interpreter, so the floor promised a
version no lane ever ran. The two settings that looked like they held it — ruff's target version
and pyright's Python version — read source syntax and not a runtime, and what a 3.12 user would
actually hit is a 3.13-only behaviour in pandas, pyarrow or numba, which neither can see. So the
floor rises to the pin, and both settings are deleted rather than retargeted: ruff derives its
target from `requires-python` and pyright takes its version from the environment's interpreter, so
Python is declared in exactly two places and declaring it a third time is how the language level
and the runtime drift apart. What it costs: a 3.12 install now fails where it used to succeed, and
this reverses a position `AGENTS.md` stated deliberately. A 3.12 test lane lost on evidence rather
than principle — nothing is known to be pinned to 3.12, and keeping the promise would mean every
bioconda pin resolving on a second interpreter — and is still the answer the day one turns up.
