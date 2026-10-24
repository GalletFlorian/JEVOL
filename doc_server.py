"""ASGI app to serve rendered markdown documentation."""

import os

import markdown

DOC_DIR = os.path.join(os.path.dirname(__file__), "doc")

PAGE_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title} - JEVOL</title>
<style>
  body {{ font-family: system-ui, sans-serif; max-width: 52rem;
         margin: 2rem auto; padding: 0 1.5rem; line-height: 1.7;
         color: #1a1a1a; }}
  h1 {{ border-bottom: 2px solid #333; padding-bottom: 0.4rem; }}
  h2 {{ margin-top: 2rem; }}
  code {{ background: #f0f0f0; padding: 0.15em 0.35em; border-radius: 3px;
          font-size: 0.92em; }}
  pre {{ background: #f4f4f4; padding: 1rem; border-radius: 6px;
         overflow-x: auto; }}
  pre code {{ background: none; padding: 0; }}
  table {{ border-collapse: collapse; margin: 1rem 0; width: 100%; }}
  th, td {{ border: 1px solid #ccc; padding: 0.5rem 0.8rem; text-align: left; }}
  th {{ background: #f0f0f0; }}
  nav {{ margin-bottom: 2rem; font-size: 0.9em; }}
  nav a {{ color: #0066cc; }}
</style>
</head>
<body>
<nav><a href="/">&larr; Index</a></nav>
{body}
</body>
</html>
"""

INDEX_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>JEVOL Documentation</title>
<style>
  body {{ font-family: system-ui, sans-serif; max-width: 600px;
         margin: 4rem auto; padding: 0 1rem; line-height: 1.6; }}
  h1 {{ border-bottom: 2px solid #333; padding-bottom: 0.5rem; }}
  ul {{ list-style: none; padding: 0; }}
  li {{ margin: 0.5rem 0; }}
  a {{ color: #0066cc; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
</style>
</head>
<body>
<h1>JEVOL Documentation</h1>
<ul>
{links}
</ul>
</body>
</html>
"""


def _find_doc_dir():
    """Find doc/ directory relative to this file or CWD."""
    candidates = [
        os.path.join(os.path.dirname(__file__), "doc"),
        os.path.join(os.getcwd(), "doc"),
    ]
    for p in candidates:
        if os.path.isdir(p):
            return os.path.abspath(p)
    raise FileNotFoundError("doc/ directory not found")


DOC_DIR = _find_doc_dir()


async def app(scope, receive, send):
    if scope["type"] != "http":
        return

    path = scope["path"].lstrip("/")

    if path == "" or path == "index.html":
        await _serve_index(send)
    elif path.endswith(".md"):
        await _serve_markdown(path, send)
    else:
        await _send_response(send, 404, b"Not found", "text/plain")


async def _serve_index(send):
    links = []
    for name in sorted(os.listdir(DOC_DIR)):
        if name.endswith(".md") and not name.startswith("_"):
            title = name.replace(".md", "").replace("_", " ").title()
            links.append(f'<li><a href="/{name}">{title}</a></li>')
    body = INDEX_TEMPLATE.format(links="\n".join(links))
    await _send_response(send, 200, body.encode(), "text/html")


async def _serve_markdown(path, send):
    filepath = os.path.join(DOC_DIR, path)
    if not os.path.isfile(filepath):
        await _send_response(send, 404, b"Not found", "text/plain")
        return

    with open(filepath, encoding="utf-8") as f:
        md_text = f.read()

    title = path.replace(".md", "").replace("_", " ").title()
    html_body = markdown.markdown(
        md_text,
        extensions=["tables", "fenced_code", "codehilite"],
    )
    page = PAGE_TEMPLATE.format(title=title, body=html_body)
    await _send_response(send, 200, page.encode(), "text/html")


async def _send_response(send, status, body, content_type):
    await send(
        {
            "type": "http.response.start",
            "status": status,
            "headers": [
                [b"content-type", content_type.encode()],
                [b"content-length", str(len(body)).encode()],
            ],
        }
    )
    await send(
        {
            "type": "http.response.body",
            "body": body,
        }
    )
