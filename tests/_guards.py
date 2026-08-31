"""The two guards every test runs behind: no network, and a **Data dir** of its own.

Both are autouse and neither can be opted out of, which is what makes "no test ever
reaches the network" and "no test ever writes into the lab's reference data" guarantees
rather than habits. They live here rather than in ``tests/conftest.py`` because a
conftest's fixtures reach the directory it sits in and nothing above it, and the
docstring examples are collected from ``src/`` — outside that tree. ``conftest.py`` at
the repository root loads this module as a plugin, which is what gives the pair the reach
the promise claims. Nothing else moved: the rest of the suite's fixtures are shared by
tests and only by tests, and stay where they were.

``install_network_guard`` is separate from the fixture that applies it so that work done
*outside* a test — a session-scoped fixture preparing a genome — can be held to the same
promise instead of being the one hole in it.
"""

from __future__ import annotations

import socket
from pathlib import Path
from typing import Any, NoReturn

import pytest
import requests

from genome.store.data_dir import LIULAB_DATA_ENV


class NetworkAccessError(RuntimeError):
    """Raised when a test attempts to reach the network."""


#: Address families that leave this machine. ``AF_UNIX`` is local IPC, never the
#: network, so the guard delegates to the real call for it.
_NETWORK_FAMILIES = frozenset({socket.AF_INET, socket.AF_INET6})

#: What to do instead — carried by every blocked call, since an error that only says
#: "no network" leaves the reader to work out how to write the test offline.
_OFFLINE_HELP = (
    "No test may reach the network. Serve a download offline with the `fake_fetch` "
    "fixture, which replaces the package's one fetch step "
    "(genome.store.fetch.fetch_url) with a copy from tests/data. A code path that "
    "also validates an assembly name at UCSC needs `requests.head` stubbed as well — "
    "see `head_recorder` in tests/assembly/test_download.py."
)


def _blocked(call: str, target: str) -> NetworkAccessError:
    """Return the error a blocked call raises: what was attempted, where, and the fix."""
    return NetworkAccessError(f"blocked network call: {call} {target}\n\n{_OFFLINE_HELP}")


def _address_text(address: Any) -> str:
    """Render a socket address as ``host:port`` when it has that shape."""
    if isinstance(address, tuple) and len(address) >= 2:
        return f"{address[0]}:{address[1]}"
    return str(address)


def install_network_guard(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make reaching the network raise, for as long as ``monkeypatch`` is in force.

    The guard itself, separated from the fixture that applies it to every test, so that
    work done *outside* a test — a session-scoped fixture preparing a genome — can be
    held to the same promise instead of being the one hole in it.
    """
    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex

    def blocked_request(
        self: requests.Session, method: Any, url: Any, *args: Any, **kwargs: Any
    ) -> NoReturn:
        raise _blocked(f"requests {str(method).upper()}", str(url))

    def blocked_connect(sock: socket.socket, address: Any) -> None:
        if sock.family in _NETWORK_FAMILIES:
            raise _blocked("socket connect", _address_text(address))
        real_connect(sock, address)

    def blocked_connect_ex(sock: socket.socket, address: Any) -> int:
        if sock.family in _NETWORK_FAMILIES:
            raise _blocked("socket connect_ex", _address_text(address))
        return real_connect_ex(sock, address)

    monkeypatch.setattr(requests.sessions.Session, "request", blocked_request)
    monkeypatch.setattr(socket.socket, "connect", blocked_connect)
    monkeypatch.setattr(socket.socket, "connect_ex", blocked_connect_ex)


@pytest.fixture(autouse=True)
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make reaching the network a failure — in every test, without being asked for.

    "No test ever touches the network" is a guarantee only if nothing can slip past
    it, so this is autouse and :func:`install_network_guard` cuts in two places:

    - ``requests.Session.request``, which every ``requests`` call funnels through: the
      ``requests.head`` of the UCSC assembly-name check, and the ``requests.get`` inside
      pooch's HTTP downloader — so a download a test forgot to replace trips here,
      naming the URL it was about to pull.
    - ``socket.socket.connect``/``connect_ex`` for the internet address families, as the
      backstop under everything that is not ``requests``: pooch's ftp and sftp
      transports, urllib, a dependency phoning home. ``AF_UNIX`` is local IPC rather
      than the network, and is left alone.

    ``pooch.retrieve`` is deliberately **not** blocked, though it is the package's only
    transport: it serves a file already sitting at the destination without any network
    call, which two tests in test_download exercise for real. Cutting at the transport
    underneath keeps those honest and still catches every download pooch would actually
    perform.

    Nothing opts out. A test that needs a download stands one in with ``fake_fetch``; a
    test that needs the UCSC name check stubs ``requests.head`` itself. Both are
    installed after this one — an autouse fixture is set up first — so they win for the
    test that asked, and the guard is back for the one that did not. A docstring example
    asks for nothing at all, which is why this reaches those too.
    """
    install_network_guard(monkeypatch)


@pytest.fixture(autouse=True)
def liulab_data(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Point the **Data dir** at this test's own directory — in every test, unasked.

    The one place the root is set. An explicit ``cache_dir`` is honoured everywhere, but
    a test exercising the *default* layout still has to point the root somewhere, and
    every such test pointing it somewhere itself is how forty-odd copies of one line
    happened. It is autouse for the same reason :func:`no_network` is: a test that forgot
    would write into the lab's real reference data, and "no test ever touches it" is a
    guarantee only if nothing can slip past.

    Request it to name the root — ``liulab_data / "genome" / "hg38"`` is where that
    assembly lands — and re-point it with ``monkeypatch.setenv`` in the few tests that
    are *about* the root: one that wants it unset, empty, or somewhere else again.

    A docstring example gets one too, and gets it for the second reason rather than the
    first: no example depends on the real root today, and this is what stops the next one
    from quietly depending on one machine's data.
    """
    monkeypatch.setenv(LIULAB_DATA_ENV, str(tmp_path))
    return tmp_path
