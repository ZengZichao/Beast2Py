#!/usr/bin/env python3
"""Fail while placeholders remain, rather than after the manuscript is submitted.

Three identifiers cannot be produced by editing: the Zenodo archive DOI (minted from a git
tag), the Dryad dataset DOI and its reviewer URL, and the TreeBASE accession.  Systematic
Biology treats an unresolved data-availability statement as an incomplete submission, and a
placeholder that survives into the deposited version is the kind of error that surfaces after
publication.  This script reads the delivered DOCX and PDF back and reports every remaining
placeholder with its location, so the check runs before upload rather than in the reviewer's
hands.

    python scripts/check_submission_placeholders.py ../Beast2Py-手稿文档
"""

from __future__ import annotations

import argparse
import re
import sys
import zipfile
from pathlib import Path

# Each pattern is a placeholder the authors must replace with a real handle.
PLACEHOLDERS = {
    "Zenodo archive DOI": re.compile(
        r"\[Zenodo[^\]]*insert[^\]]*\]", re.I),
    "Dryad dataset DOI": re.compile(r"\[Dryad (?:dataset|DOI)[^\]]*insert[^\]]*\]", re.I),
    "Dryad reviewer URL": re.compile(r"\[Dryad reviewer[^\]]*\]", re.I),
    "TreeBASE accession": re.compile(r"\[TreeBASE[^\]]*insert[^\]]*\]", re.I),
    "generic bracketed gap": re.compile(
        r"\[(?:dataset|DOI|URL|accession|ID)[^\]]*insert[^\]]*\]", re.I),
    "unresolved TODO": re.compile(r"\b(?:TODO|TBD|XXX+)\b"),
}


def paragraphs_of(docx_path: Path):
    """Yield (index, text) for each body paragraph, without needing python-docx."""
    with zipfile.ZipFile(docx_path) as z:
        xml = z.read("word/document.xml").decode("utf-8")
    for n, para in enumerate(re.findall(r"<w:p[ >].*?</w:p>", xml, re.S), 1):
        text = "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", para, re.S))
        text = (text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">"))
        if text.strip():
            yield n, text


def pdf_pages_of(pdf_path: Path):
    import fitz

    doc = fitz.open(str(pdf_path))
    for n, page in enumerate(doc, 1):
        yield n, page.get_text()
    doc.close()


def scan(folder: Path) -> list:
    problems = []
    for docx in sorted(folder.glob("*.docx")):
        for n, text in paragraphs_of(docx):
            for label, pat in PLACEHOLDERS.items():
                for m in pat.finditer(text):
                    problems.append("%s  para %d  [%s]  …%s…"
                                    % (docx.name, n, label,
                                       text[max(0, m.start() - 40):m.end() + 40].strip()))
    for pdf in sorted(folder.glob("*.pdf")):
        try:
            pages = list(pdf_pages_of(pdf))
        except ImportError:
            print("  (skipping %s: pymupdf not installed)" % pdf.name)
            continue
        for n, text in pages:
            for label, pat in PLACEHOLDERS.items():
                if pat.search(text):
                    problems.append("%s  page %d  [%s]" % (pdf.name, n, label))
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folders", nargs="*", type=Path,
                    help="folders holding the manuscript DOCX/PDF")
    ap.add_argument("--expect", type=int, default=None,
                    help="assert exactly this many open placeholders and fail otherwise")
    args = ap.parse_args(argv)
    if not args.folders:
        ap.error("name at least one folder")

    problems = []
    for folder in args.folders:
        problems += scan(folder)
    for line in problems:
        print("  - %s" % line)
    print("%d unresolved identifier(s) across %d file(s)"
          % (len(problems), sum(len(list(f.glob('*.docx'))) + len(list(f.glob('*.pdf')))
                                for f in args.folders)))
    if args.expect is not None and len(problems) != args.expect:
        print("expected %d, found %d" % (args.expect, len(problems)))
        return 1
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
