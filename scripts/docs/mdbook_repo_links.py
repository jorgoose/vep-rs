#!/usr/bin/env python3
"""mdBook preprocessor that makes the repository's markdown read correctly as a site.

The book's source is `docs/`, but those pages link outside it (`../README.md`,
`../scripts/README.md`, the provenance JSON) and the site's home page is the
repository README, whose links are relative to the repository root. mdBook
copies only `docs/` into the site, so those links would 404. This preprocessor:

- replaces the content of the `index.md` chapter with `README.md`, so the
  README is the single source of the home page;
- resolves every relative markdown link against the file it was written in;
- keeps a link to a page in the book relative (mdBook turns `.md` into `.html`),
  points `README.md` at the home page, and sends any other repository path to
  its GitHub view on `main`.

Links inside fenced code blocks are left alone. The GitHub base URL comes from
`output.html.git-repository-url` in `book.toml`.

mdBook runs it (`[preprocessor.repo-links]` in `book.toml`): `supports <renderer>`
exits 0, and otherwise it reads `[context, book]` JSON on stdin and writes the
book back on stdout. Stdlib only.
"""

import json
import os
import re
import sys

LINK = re.compile(r"(\]\()([^)\s]+)((?:\s+\"[^\"]*\")?\))")
FENCE = re.compile(r"^\s*(```|~~~)")
SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
INDEX = "index.md"


def rewrite_target(target, base_dir, chapter_dir, src, root, chapters, blob_url):
    """Return the site-correct form of one link target.

    base_dir is the repository-relative directory the link was written in;
    chapter_dir is the chapter's directory relative to the book source.
    """
    if target.startswith("#") or SCHEME.match(target):
        return target
    path, sep, frag = target.partition("#")
    if not path:
        return target
    resolved = os.path.normpath(os.path.join(base_dir, path)).replace(os.sep, "/")
    if resolved.startswith("../"):
        return target
    in_book = None
    if resolved == "README.md":
        in_book = INDEX
    elif resolved.startswith(src + "/") and resolved[len(src) + 1 :] in chapters:
        in_book = resolved[len(src) + 1 :]
    if in_book is not None:
        rel = os.path.relpath(in_book, chapter_dir or ".").replace(os.sep, "/")
        return rel + sep + frag
    kind = "tree" if os.path.isdir(os.path.join(root, resolved)) else "blob"
    return f"{blob_url}/{kind}/main/{resolved}" + sep + frag


def rewrite(content, **kw):
    out, in_fence = [], False
    for line in content.splitlines(keepends=True):
        if FENCE.match(line):
            in_fence = not in_fence
        elif not in_fence:
            line = LINK.sub(lambda m: m[1] + rewrite_target(m[2], **kw) + m[3], line)
        out.append(line)
    return "".join(out)


def walk(items):
    for item in items:
        chapter = item.get("Chapter") if isinstance(item, dict) else None
        if chapter:
            yield chapter
            yield from walk(chapter.get("sub_items", []))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "supports":
        sys.exit(0)
    context, book = json.load(sys.stdin)
    root = context["root"]
    config = context["config"]
    src = os.path.normpath(config.get("book", {}).get("src", "src")).replace(os.sep, "/")
    repo_url = config.get("output", {}).get("html", {}).get("git-repository-url")
    if not repo_url:
        sys.exit("mdbook_repo_links: output.html.git-repository-url is not set in book.toml")
    blob_url = repo_url.rstrip("/")

    chapters = {c["source_path"] for c in walk(book["items"]) if c.get("source_path")}
    for chapter in walk(book["items"]):
        source = chapter.get("source_path")
        if not source:
            continue
        chapter_dir = os.path.dirname(source)
        if source == INDEX:
            with open(os.path.join(root, "README.md"), encoding="utf-8") as f:
                chapter["content"] = f.read()
            base_dir = ""
        else:
            base_dir = os.path.join(src, chapter_dir) if chapter_dir else src
        chapter["content"] = rewrite(
            chapter["content"],
            base_dir=base_dir,
            chapter_dir=chapter_dir,
            src=src,
            root=root,
            chapters=chapters,
            blob_url=blob_url,
        )
    json.dump(book, sys.stdout)


if __name__ == "__main__":
    main()
