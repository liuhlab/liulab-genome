"""The package's one fetch step — every byte it downloads comes through here.

A module of its own because two things fetch: :mod:`genome.io.download`, which brings an
assembly's FASTA over the network, and :mod:`genome.io.annotation.registration`, which
brings an annotation's
GTF. While the step lived under the downloader the second of those was an import back
into it, and that was the last edge of a cycle — ``download`` reaches ``chimera`` to
build one, ``chimera`` reaches ``gtf`` to merge the annotations, and ``gtf`` reached
``download`` for this one call.

**Deliberately a leaf**: it imports nothing else from this package, which is what lets
both callers have it without either having the other. Nothing here knows what an
assembly, an annotation or a **Completion marker** is — it is given a URL and a
directory, and it returns a path.

Examples
--------
>>> from genome.io import fetch
>>> fetch.fetch_url                                      # doctest: +ELLIPSIS
<function fetch_url at ...>
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pooch

# A pooch post-processor: called with (fname, action, pooch_instance) and
# returns the path (or paths) to use as the result of the download.
Processor = Callable[..., object]


def fetch_url(
    url: str,
    dest_dir: Path,
    *,
    known_hash: str | None = None,
    fname: str | None = None,
    processor: Processor | None = None,
    progressbar: bool = True,
) -> Path:
    """Download ``url`` into ``dest_dir`` and return the local path.

    The package's one fetch step: every download it performs goes through here, and
    this is the only call site of :func:`pooch.retrieve`. pooch picks its transport
    from the URL scheme, so ``http``, ``https``, ``ftp`` and ``sftp`` all work — the
    last additionally needs ``paramiko`` installed, which pooch will say for itself.
    A file already sitting at the destination is reused (verified against
    ``known_hash`` when one is given) and no network call is made.

    Reach this function through the module, never by importing the name: write
    ``from genome.io import fetch`` and call ``fetch.fetch_url(...)``, so that a single
    ``monkeypatch.setattr(fetch, "fetch_url", ...)`` takes every download in the package
    offline. Both callers — the assembly downloader and the annotation registration —
    spell it that way for that reason, and the suite's offline guard rests on it: a
    caller that imported the name instead would hold a reference the rebinding never
    reaches, and would keep downloading while the tests believed themselves offline.

    Parameters
    ----------
    url : str
        The file URL to download, including its scheme.
    dest_dir : pathlib.Path
        Directory the file is written into. Created when a download actually happens.
    known_hash : str, optional
        Expected hash as ``"<algorithm>:<hexdigest>"`` (e.g. ``"md5:8f3c..."``), or a
        bare hex digest for sha256. If ``None``, verification is skipped and pooch logs
        the computed hash so you can pin it next time.
    fname : str, optional
        Local file name to save as. Defaults to a hash-prefixed unique name pooch
        derives from ``url``.
    processor : callable, optional
        A pooch post-processor applied after the download, such as
        :class:`pooch.Decompress` or :class:`pooch.Untar`. Its return value becomes the
        path returned here.
    progressbar : bool, default True
        Show a textual download progress bar (requires ``tqdm``).

    Returns
    -------
    pathlib.Path
        Absolute path to the downloaded (and, if ``processor`` was given, processed)
        file.

    Raises
    ------
    requests.exceptions.HTTPError
        If an http(s) download fails (e.g. the URL 404s).
    ValueError
        If ``known_hash`` is given and the file does not match it, or if no downloader
        exists for the URL's scheme.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io import fetch
    >>> fetch.fetch_url(                                     # doctest: +SKIP
    ...     "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
    ...     Path("/scratch/liulab/genome/sacCer3"),
    ...     fname="sacCer3.fa.gz",
    ... )
    PosixPath('/scratch/liulab/genome/sacCer3/sacCer3.fa.gz')
    """
    result = pooch.retrieve(
        url=url,
        known_hash=known_hash,
        fname=fname,
        path=dest_dir,
        processor=processor,
        progressbar=progressbar,
    )
    return Path(result)
