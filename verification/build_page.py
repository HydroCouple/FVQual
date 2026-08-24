#!/usr/bin/env python3
"""Inline verification/data.json into verification/page.html -> index.html.

The page is a single self-contained file so it can be opened from disk,
committed, or published without a server or any external asset.
"""
import json, pathlib, sys

here = pathlib.Path(__file__).parent
data = json.loads((here / "data.json").read_text())

# Trim float precision: animation frames do not need 6 significant digits,
# and this roughly halves the page.
def trim(o):
    if isinstance(o, float):
        return round(o, 5)
    if isinstance(o, list):
        return [trim(v) for v in o]
    if isinstance(o, dict):
        return {k: trim(v) for k, v in o.items()}
    return o

payload = json.dumps(trim(data), separators=(",", ":"))
# The payload sits in a <script type="application/json"> block, so the only
# sequence that could break out is a literal closing script tag.
payload = payload.replace("</", "<\\/")

html = (here / "page.html").read_text().replace("__DATA__", payload)
out = here / "index.html"
out.write_text(html)

print(f"wrote {out} ({out.stat().st_size/1e6:.2f} MB)")
