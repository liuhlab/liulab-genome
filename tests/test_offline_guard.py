"""The suite-wide offline guard: reaching the network fails the test that tried.

``no_network`` in conftest is autouse, so every test in this suite already runs behind
it; these are what prove it fires rather than merely being installed. Each deliberately
attempts a call that would leave this machine and asserts it raises instead. None of
them requests the guard, which is the point — it is on without being asked for.
"""

from __future__ import annotations

import socket
from pathlib import Path

import pytest
import requests

from genome.io import fetch as fetch_mod
from genome.io.download import UCSCGenomeDownloader

from .conftest import NetworkAccessError

#: A real URL, so a guard that failed to fire would be visible as a live download.
_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"


@pytest.mark.parametrize("verb", ["head", "get", "post"])
def test_a_requests_call_is_blocked_and_says_what_to_do_instead(verb: str) -> None:
    with pytest.raises(NetworkAccessError) as excinfo:
        getattr(requests, verb)(_URL)

    message = str(excinfo.value)
    assert verb.upper() in message  # which call was attempted
    assert _URL in message  # and to what
    assert "fake_fetch" in message  # and what serves it offline instead


def test_the_packages_one_fetch_step_is_blocked(tmp_path: Path) -> None:
    # Nothing at the destination, so pooch has to download it — through requests, into
    # the guard. This is what a test that forgot ``fake_fetch`` would hit.
    with pytest.raises(NetworkAccessError, match=r"hg38\.fa\.gz"):
        fetch_mod.fetch_url(_URL, tmp_path, fname="hg38.fa.gz", progressbar=False)


def test_the_ucsc_assembly_name_check_is_blocked(tmp_path: Path) -> None:
    # The gap this guard closes: an assembly whose metadata pins no source URL has its
    # name validated with a live HEAD at UCSC, and nothing patched that repo-wide.
    downloader = UCSCGenomeDownloader("hg38", tmp_path)

    with pytest.raises(NetworkAccessError, match="requests HEAD"):
        downloader.validate_assembly()


@pytest.mark.parametrize("method", ["connect", "connect_ex"])
def test_a_tcp_connection_is_blocked(method: str) -> None:
    # The backstop under everything that is not requests — pooch's ftp/sftp transports,
    # urllib, a dependency phoning home. 192.0.2.1 is TEST-NET-1 and routes nowhere
    # anyway; the guard raises before the address is so much as resolved.
    with (
        socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock,
        pytest.raises(NetworkAccessError, match=r"192\.0\.2\.1:80"),
    ):
        getattr(sock, method)(("192.0.2.1", 80))


def test_a_local_unix_socket_is_left_alone(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # The one thing the socket cut deliberately does not block: AF_UNIX is local IPC,
    # not the network. A relative name keeps the path under the AF_UNIX length limit.
    monkeypatch.chdir(tmp_path)
    with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as server:
        server.settimeout(5)
        server.bind("s")
        server.listen(1)
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
            client.settimeout(5)
            client.connect("s")
            connection, _ = server.accept()
            with connection:
                client.sendall(b"ping")
                assert connection.recv(4) == b"ping"
