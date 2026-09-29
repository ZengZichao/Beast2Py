#!/usr/bin/env python3
"""Gate the deliverables that no Python script redraws.

Every figure in this submission except Figures 1, 2 and S4 is redrawn from a
``.py`` source that carries its own compliance checks, so a bad export fails at
draw time.  The three draw.io figures have no such hook: the ``.drawio`` file is
the source and the ``.pdf``/``.png`` pair is hand-exported from a GUI, which is
how Figure 1 shipped as a PDF 19 hours older than the source it was exported
from.  This script reads the delivered artefacts back and applies the same
thresholds the Python figures apply to themselves.

Run it after every draw.io export::

    python scripts/check_figure_export.py ../Beast2Py-手稿正图 ../Beast2Py-手稿附图

Optional dependencies: pymupdf, pillow (``pip install -e .[figures]``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# The same floors the figure scripts gate themselves against, restated here
# rather than imported: those files are submission sources, not a package.
TEXT_MIN_PT = 7.92
HAIRLINE_MIN_PT = 0.71
PPI_MIN = 600.0
ALLOWED_FONTS = {"Helvetica", "Helvetica-Bold", "Helvetica-Oblique"}


def _fail(problems, msg):
    problems.append(msg)


def check_pdf(pdf, problems):
    """Return (page_width_in, smallest text pt) for one exported PDF."""
    import fitz

    doc = fitz.open(str(pdf))
    if len(doc) != 1:
        _fail(problems, "%s: %d pages, a figure is one page" % (pdf.name, len(doc)))
    page = doc[0]
    width_in = page.rect.width / 72.0

    sizes, fonts = [], set()
    for block in page.get_text("dict")["blocks"]:
        if block.get("type") == 1:
            _fail(problems, "%s: page %d carries a raster image object"
                            % (pdf.name, block.get("number", 0)))
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                if not span["text"].strip():
                    continue
                sizes.append(span["size"])
                fonts.add(span["font"])
    for name in sorted(fonts - ALLOWED_FONTS):
        _fail(problems, "%s: font %r is outside the Helvetica family"
                        % (pdf.name, name))
    if "Type3" in fonts:
        _fail(problems, "%s: Type3 fonts are not embeddable for print" % pdf.name)

    thin = HAIRLINE_MIN_PT
    seen_rule = False
    for drawing in page.get_drawings():
        width = drawing.get("width") or 0.0
        if width and drawing.get("color") is not None:
            seen_rule = True
            thin = min(thin, width)
    if not seen_rule:
        _fail(problems, "%s: no stroked path found, so nothing was drawn" % pdf.name)
    if thin < HAIRLINE_MIN_PT:
        _fail(problems, "%s: thinnest rule %.3f pt, under %.2f pt"
                        % (pdf.name, thin, HAIRLINE_MIN_PT))

    smallest = min(sizes) if sizes else 0.0
    if sizes and smallest < TEXT_MIN_PT - 1e-9:
        _fail(problems, "%s: text at %.2f pt, under the %.2f pt floor"
                        % (pdf.name, smallest, TEXT_MIN_PT))
    doc.close()
    return width_in, smallest, thin


def check_png(png, page_width_in, problems):
    from PIL import Image

    im = Image.open(str(png))
    ppi = im.size[0] / page_width_in if page_width_in else 0.0
    if ppi < PPI_MIN:
        _fail(problems, "%s: %.0f px over a %.3f in page is %.0f ppi, under %g"
                        % (png.name, im.size[0], page_width_in, ppi, PPI_MIN))
    return im.size[0], ppi


def check_figure(drawio: Path, problems: list) -> str:
    for suffix in (".pdf", ".png"):
        out = drawio.with_suffix(suffix)
        if not out.exists():
            _fail(problems, "%s: no delivered %s beside this source"
                            % (drawio.name, suffix.lstrip(".")))
            return ""
        if out.stat().st_mtime < drawio.stat().st_mtime:
            _fail(problems, "%s is older than %s: re-export it"
                            % (out.name, drawio.name))
    width_in, smallest, thin = check_pdf(drawio.with_suffix(".pdf"), problems)
    px, ppi = check_png(drawio.with_suffix(".png"), width_in, problems)
    return ("%s: page %.3f in, smallest text %.2f pt, thinnest rule %.3f pt, "
            "png %d px = %.0f ppi" % (drawio.stem, width_in, smallest, thin, px, ppi))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folders", nargs="*", type=Path,
                    help="folders to scan for .drawio sources")
    args = ap.parse_args(argv)
    if not args.folders:
        ap.error("name at least one folder holding .drawio sources")

    sources = sorted(p for folder in args.folders for p in folder.glob("*.drawio"))
    if not sources:
        print("no .drawio sources in %s" % ", ".join(str(f) for f in args.folders))
        return 0

    problems: list = []
    for src in sources:
        line = check_figure(src, problems)
        if line:
            print(line)
    for problem in problems:
        print("  !! %s" % problem)
    if problems:
        print("%d problem(s) over %d draw.io figure(s)"
              % (len(problems), len(sources)))
        return 1
    print("ok: all %d draw.io deliverables pass the same floors the script-drawn "
          "figures apply to themselves" % len(sources))
    return 0


if __name__ == "__main__":
    sys.exit(main())
