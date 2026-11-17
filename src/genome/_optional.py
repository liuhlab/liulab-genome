"""Libraries the package uses but does not declare, imported where they are used.

``MOODS`` and ``gffutils`` are not in ``[project] dependencies``: PyPI carries MOODS only as
``MOODS-python``, an sdist that needs SWIG and a C++ compiler, so declaring it would turn a
working ``pip install`` into a failing build. ``import genome`` therefore has to work without
them, and a call that needs one fails at call time, saying how to get it.

The import itself stays a plain ``import`` statement inside the function that uses the
library, so the one-entrance tests and the type checker still see it. What lives here is the
error that follows a failed one.
"""

from __future__ import annotations

#: How to get MOODS. Pip has only ``MOODS-python``, which compiles C++ on install.
MOODS_INSTALL = "Install `moods` from bioconda; pip's `MOODS-python` needs a C++ compiler."

#: How to get gffutils.
GFFUTILS_INSTALL = "Install `gffutils`."


def missing_library(name: str, install: str) -> ModuleNotFoundError:
    """Build the error for a call that needs a library that is not installed.

    Raise it from the ``ImportError`` that reported the library missing, so the original
    stays attached as the cause.

    Parameters
    ----------
    name : str
        The library, as its import name.
    install : str
        What the caller should do to get it, as one sentence.

    Returns
    -------
    ModuleNotFoundError
        Naming the library and ``install``.

    Examples
    --------
    >>> str(missing_library("gffutils", GFFUTILS_INSTALL))
    'gffutils is not installed, and this call needs it. Install `gffutils`.'
    """
    return ModuleNotFoundError(
        f"{name} is not installed, and this call needs it. {install}", name=name
    )
