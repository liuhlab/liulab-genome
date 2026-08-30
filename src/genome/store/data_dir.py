"""Where the lab's reference data begins, and the roots filed directly under it.

The **Data dir** belongs to no context — the assembly tree is most of it, a **Prepared
set** is filed beside that tree, and every context reads one or the other — so the
environment variable that names it is read here and nowhere else. What each root then
*is* belongs to whoever fills it: :func:`~genome.assembly.registration.assembly_data_dir`
is the Assembly context's, and the three prepared-set roots are the Motif, Xref and
Orthology contexts' own, each declared beside the code that reads it.

Nothing here creates a directory. A path is an answer to *where would this go*, and the
write that needs it is what brings it into existence.

Examples
--------
>>> import os
>>> from genome.store.data_dir import liulab_data_dir, prepared_data_dir
>>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
>>> liulab_data_dir()
PosixPath('/scratch/liulab')
>>> prepared_data_dir("xref")
PosixPath('/scratch/liulab/xref')
>>> del os.environ["LIULAB_DATA"]
"""

from __future__ import annotations

import os
from pathlib import Path

#: Environment variable naming the lab data root directory.
LIULAB_DATA_ENV = "LIULAB_DATA"

#: Well-known lab data roots, tried in order when ``LIULAB_DATA`` is unset.
DEFAULT_LIULAB_DATA_PATHS = [
    "/share/lhqlab/liulab_data",
    "/large_storage/zhoulab/hanliu/liulab_data",
]


def liulab_data_dir() -> Path:
    """Return the root directory for lab reference data.

    The location is read from the ``LIULAB_DATA`` environment variable. When that
    is unset (or empty), each entry in :data:`DEFAULT_LIULAB_DATA_PATHS` is checked
    in order and the first that exists is used as the root. If none exist, it falls
    back to ``~/liulab_data``. The path is expanded (``~`` resolved) but **not**
    created here — callers create the specific subdirectory they need on first write.

    Returns
    -------
    pathlib.Path
        The resolved lab data root.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> liulab_data_dir()
    PosixPath('/scratch/liulab')
    >>> del os.environ["LIULAB_DATA"]
    """
    env = os.environ.get(LIULAB_DATA_ENV)
    if env:
        return Path(env).expanduser()
    for candidate in DEFAULT_LIULAB_DATA_PATHS:
        path = Path(candidate).expanduser()
        if path.exists():
            return path
    return (Path.home() / "liulab_data").expanduser()


def prepared_data_dir(subdir: str) -> Path:
    """Return one **Prepared set** root under the **Data dir**, created by nobody.

    A prepared set is a *sibling* of the assembly tree rather than a tenant of it —
    ``motif/``, ``xref/``, ``homology/`` — because none of the three belongs to an
    **Assembly**. Which name a context uses, and what lives under it, is that context's
    own business and is declared beside the code that reads it; this is the one line all
    three share.

    Parameters
    ----------
    subdir : str
        The root's directory name, as its own context spells it.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/<subdir>``. Nothing is created by asking.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> prepared_data_dir("motif")
    PosixPath('/scratch/liulab/motif')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / subdir
