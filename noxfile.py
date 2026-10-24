"""Nox sessions for JEVOL."""

import nox

nox.options.sessions = ["tests"]


@nox.session(python=["3.12", "3.13"])
def dev(session: nox.Session) -> None:
    """Install development dependencies."""
    session.run("uv", "sync", "--extra", "dev", external=True)


@nox.session(python=["3.12", "3.13"])
def tests(session: nox.Session) -> None:
    """Run the test suite."""
    session.run("uv", "sync", "--extra", "dev", external=True)
    session.run("uv", "run", "pytest", "-v", *session.posargs, external=True)


@nox.session(python="3.13")
def coverage(session: nox.Session) -> None:
    """Run tests with coverage."""
    session.run("uv", "sync", "--extra", "dev", external=True)
    session.run(
        "uv",
        "run",
        "pytest",
        "--cov=Python",
        "--cov-report=term-missing",
        "--cov-report=xml",
        "-v",
        *session.posargs,
        external=True,
    )


@nox.session(python="3.13")
def lint(session: nox.Session) -> None:
    """Run linting with ruff."""
    session.run("uv", "sync", "--extra", "dev", external=True)
    session.run("uv", "run", "ruff", "check", "Python", "noxfile.py", "doc_server.py", *session.posargs, external=True)


@nox.session(python="3.13")
def format(session: nox.Session) -> None:
    """Format code with ruff."""
    session.run("uv", "sync", "--extra", "dev", external=True)
    session.run("uv", "run", "ruff", "format", "Python", "noxfile.py", "doc_server.py", *session.posargs, external=True)


@nox.session(python="3.13")
def serve_doc(session: nox.Session) -> None:
    """Serve rendered markdown docs at http://localhost:8000."""
    session.run("uv", "sync", "--extra", "doc", external=True)
    session.run("uv", "run", "uvicorn", "doc_server:app", "--host", "0.0.0.0", "--port", "8000", external=True)
