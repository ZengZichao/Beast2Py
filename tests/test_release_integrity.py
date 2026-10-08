"""Release-integrity guards.

Each test locks a defect that must not silently return:

* the birthday-collision bound quoted to users,
* comment text that lost a character to an earlier edit pass,
* the sidecar/XML pair that lets a reader verify an analysis without running
  anything,
* order-stability of the fingerprint serialiser,
* what the built wheel and sdist must contain, checked from the declarations,
* the tutorial steps that pointed two commands at one output name, and the
  documented commands a reader is told to run.
"""

import json
import math
import pathlib
import re

import pytest

import beast2py
from beast2py.reproducibility import FingerprintGenerator

PKG = pathlib.Path(beast2py.__file__).parent
REPO = PKG.parent
EXAMPLES = REPO / "examples"


class TestCollisionBound:
    """The 48-bit identifier's collision caveat must stay arithmetically right."""

    def test_hash_carries_48_bits(self):
        assert FingerprintGenerator.collision_bits() == 48

    def test_fifty_percent_bound_is_not_the_square_root_of_the_space(self):
        bits = FingerprintGenerator.collision_bits()
        space = 2**bits
        n50 = FingerprintGenerator.collision_bound_50pct()
        # P(at least one collision) at the reported n must really be ~0.5
        p = 1 - math.exp(-n50 * (n50 - 1) / (2 * space))
        assert abs(p - 0.5) < 1e-4, p
        # sqrt(N) alone is only a 39% chance; conflating the two was the defect.
        p_sqrt = 1 - math.exp(-math.sqrt(space) * (math.sqrt(space) - 1) / (2 * space))
        assert p_sqrt < 0.45

    def test_note_quotes_the_bound_and_not_2_pow_24(self):
        note = FingerprintGenerator.collision_note()
        assert "2**24" not in note
        assert "1.98e+07" in note
        assert "change-detection token" in note

    def test_committed_sidecar_carries_the_corrected_note(self):
        sidecar = EXAMPLES / "output_basic.fingerprint.json"
        if not sidecar.exists():
            pytest.skip("example sidecar not present in this checkout")
        note = json.loads(sidecar.read_text(encoding="utf-8"))["collision_note"]
        assert "2**24" not in note
        assert "ln 2" in note


class TestSidecarPairing:
    """A sidecar must describe the XML it sits beside, or it is worse than none."""

    def test_committed_pair_is_self_consistent(self):
        import hashlib

        xml = (EXAMPLES / "output_basic.xml").read_bytes()
        side = json.loads((EXAMPLES / "output_basic.fingerprint.json").read_text())
        assert side["xml_digest"] == hashlib.sha256(xml).hexdigest()
        embedded = re.search(
            rb"Analysis Fingerprint: (B2P-[0-9a-f]{12}-\S+?) \| Data: ([0-9a-f]+)", xml
        )
        assert embedded, "XML lost its fingerprint comment"
        assert embedded.group(1).decode() == side["fingerprint"]
        assert embedded.group(2).decode() == side["data_hash"]
        assert side["config_hash"].startswith(side["fingerprint"].split("-")[1])

    def test_committed_trio_is_what_the_shipped_code_produces(self):
        """Self-consistency is not the same thing as being current.

        The XML and its sidecar can age together behind a change in the writer, and
        the reader who runs the documented command would then get something other
        than the committed example the documentation presents as validated.
        """
        import contextlib
        import io
        import tempfile

        from beast2py.api import generate_xml
        from beast2py.config import ConfigParser

        out = pathlib.Path(tempfile.mkdtemp()) / "output_basic.xml"
        with contextlib.redirect_stdout(io.StringIO()):
            generate_xml(str(EXAMPLES / "config_basic.yaml"), str(out), force=True)
        assert (
            out.read_bytes() == (EXAMPLES / "output_basic.xml").read_bytes()
        ), "the committed exemplar XML is behind the writer"

        config = ConfigParser().parse(EXAMPLES / "config_basic.yaml")
        side = json.loads((EXAMPLES / "output_basic.fingerprint.json").read_text())
        fresh = FingerprintGenerator.generate_fingerprint_dict(
            config, xml_content=out.read_text(encoding="utf-8")
        )
        volatile = {"generation_time"}
        assert {k: v for k, v in side.items() if k not in volatile} == {
            k: v for k, v in fresh.items() if k not in volatile
        }, "the committed sidecar is behind what the shipped code writes"

    def test_fingerprint_literals_in_docs_are_still_producible(self):
        """Every digest the documentation shows must be one the tool emits today.

        The XML-format page carried the real alignment digest together with a
        configuration digest no example produces any more, which is exactly the
        shape a reader would take as evidence that the fingerprint scheme is
        reproducible -- and then fail to reproduce.
        """
        import contextlib
        import io

        from beast2py.config import ConfigParser

        cfg_hashes, data_hashes = set(), set()
        for cfg in sorted(EXAMPLES.glob("config_*.yaml")):
            buf = io.StringIO()
            try:
                with contextlib.redirect_stdout(buf):
                    parsed = ConfigParser().parse(cfg)
            except Exception:  # noqa: BLE001 - a config this suite cannot parse is locked elsewhere
                continue
            cfg_hashes.add(FingerprintGenerator.config_hash(parsed))
            data_hashes.add(FingerprintGenerator.data_hash(parsed))

        stale = []
        for md in sorted(REPO.rglob("*.md")):
            if "subprojects" in str(md):
                continue
            text = md.read_text(encoding="utf-8", errors="replace")
            for digest in re.findall(r"B2P-([0-9a-f]{12})-", text):
                if digest not in cfg_hashes:
                    stale.append("%s: config digest %s" % (md.name, digest))
            for digest in re.findall(r"(?:Data: |\"data_hash\": \")([0-9a-f]{16})", text):
                if digest not in data_hashes:
                    stale.append("%s: data digest %s" % (md.name, digest))
        assert not stale, "documentation shows digests the code cannot produce: %s" % stale[:6]


class TestCanonicalSerialisation:
    """A set must not be able to make the digest depend on PYTHONHASHSEED."""

    def test_sets_serialise_in_sorted_order(self):
        a = FingerprintGenerator._canonical({"z", "a", "m"})
        b = FingerprintGenerator._canonical({"m", "z", "a"})
        assert a == b == ["a", "m", "z"]

    def test_dicts_serialise_in_sorted_order(self):
        import json

        one = json.dumps({"b": 1, "a": 2}, default=FingerprintGenerator._canonical, sort_keys=True)
        two = json.dumps({"a": 2, "b": 1}, default=FingerprintGenerator._canonical, sort_keys=True)
        assert one == two


class TestCommentIntegrity:
    """Earlier edit passes stripped characters from comments and docstrings.

    The failures were mechanical, so they are mechanically detectable -- but only
    by a guard that looks where they actually happened. The first version of this
    class missed precisely the shapes it existed to catch: it read only hash
    comments and never a docstring, it required four or more leading spaces for a
    prose line (so a module docstring line at column 0 was invisible), and it
    counted parentheses only inside comment blocks (so an unbalanced docstring was
    invisible). Those blind spots are closed here, and each new rule is pinned by a
    test_guard method that writes a broken sample and asserts the matching predicate
    fires. A guard that cannot fail is worse than no guard: it reads like a record
    that the defect class was checked.
    """

    # The package sources *and* the test tree: the same edit passes touched both,
    # and a guard that watches only the shipped package leaves the harness
    # unmonitored.
    TARGETS = sorted(p for root in (PKG, REPO / "tests") for p in root.rglob("*.py"))

    # A period preceded by a space is the signature of an edit that deleted the
    # word after it. The alternatives are what a repaired sentence legitimately
    # ends with: end of line, the opening of the next sentence, a docstring
    # terminator, or a dash. "from . import" stays unmatched because "import"
    # begins with a lower-case letter.
    ANCHOR = re.compile(r"(\w|[^\x00-\x7f]) +\.(?:$|\"| [A-Z`\"'#]| [-\u2014])")
    # A line that is nothing but a period, alone or behind a hash sign.
    EMPTY = re.compile(r"^\s*(?:#\s*)?\.\s*$")
    # A line that is nothing but a numbered anchor: the residue of a list whose
    # items were merged into prose and whose numbering was left behind.
    STRAY_ANCHOR = re.compile(r"^\s*#?\s*\((?:\d+|[a-z])\)\s*\.?\s*$")
    # An opening parenthesis followed by a space is what a word swallowed right
    # after it looks like. The escape keeps this file from tripping its own rule.
    OPEN_PAREN = re.compile(r"\x28 \S")
    # A numbered anchor is only meaningful as part of a list that starts at one
    # and does not skip.
    NUMBERED = re.compile(r"\((\d+)\)")
    # A *sentence* closing on a bare numbered anchor is the other half of the same
    # defect: the list its item belonged to was merged into prose and only the
    # numbering survived. The sequence rule above cannot see it, because a lone
    # anchor is trivially "a list that starts at one", so the shape needs its own
    # test. Built from fragments: this file is one of the files the scan covers.
    TAIL_ANCHOR = re.compile(r"\x28\d+\x29\s*(?=\.\s|\.$|\Z)")
    # Function words that cannot close a sentence: their presence at the end of a
    # docstring or comment block means an edit cut the text short, which is the
    # dominant defect class in this repository.
    DANGLING = re.compile(
        r"\b(?:the|a|an|of|and|or|is|was|for|with|from|which|to|in|by|as)\s*$",
        re.IGNORECASE,
    )

    @staticmethod
    def _lines(path):
        return path.read_text(encoding="utf-8").split("\n")

    @classmethod
    def _prose_line_numbers(cls, path):
        """Return the line numbers that hold prose: comments and docstring bodies.

        Code is deliberately out of scope. Half of the repository builds BEAST2
        XML and LaTeX out of multi-line string literals, and those carry the
        punctuation of the *generated* document, not of a sentence an edit pass
        can truncate.
        """
        wanted = set()
        for start, text in cls._docstrings(path):
            wanted.update(range(start, start + len(text.split("\n"))))
        for no, line in enumerate(cls._lines(path), 1):
            if line.strip().startswith("#"):
                wanted.add(no)
        return wanted

    @classmethod
    def _prose_lines(cls, path):
        """Yield every prose line with its number, with no indentation floor.

        The floor the first version used excluded exactly the lines that were
        damaged: a module docstring and a top-of-file comment both sit at column 0.
        """
        wanted = cls._prose_line_numbers(path)
        for no, line in enumerate(cls._lines(path), 1):
            if no in wanted and line.strip():
                yield no, line

    @staticmethod
    def _docstrings(path):
        """Yield (first_line, text) for every module, class and function docstring."""
        import ast

        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(
                node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            ):
                continue
            body = getattr(node, "body", [])
            if not body:
                continue
            first = body[0]
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                yield first.value.lineno, first.value.value

    @classmethod
    def _prose_units(cls, path):
        """Yield (line, joined_text) per docstring and per contiguous comment block.

        Sentences are judged on the joined unit, never on one wrapped line: a
        comment that wraps mid-sentence is normal, one that ends mid-sentence is
        the damage.
        """
        for start, text in cls._docstrings(path):
            joined = " ".join(line.strip() for line in text.split("\n") if line.strip())
            if joined:
                yield start, joined
        block, start = [], 0
        for no, line in enumerate(cls._lines(path), 1):
            if line.strip().startswith("#"):
                if not block:
                    start = no
                block.append(line.strip().lstrip("#").strip())
                continue
            if block:
                yield start, " ".join(block)
                block = []
        if block:
            yield start, " ".join(block)

    @classmethod
    def paren_counts(cls, text):
        """Return the open and close counts of `text`, ignoring interval notation.

        "(0, 1]" and "[1, 2)" legitimately unbalance, and this repository writes
        about parameter bounds, so bracketed numeric intervals are removed first.
        """
        cleaned = re.sub(r"[\[(]\s*[\d.,\s-]+?[\])]", "", text)
        return cleaned.count("("), cleaned.count(")")

    @classmethod
    def stray_numbered_anchors(cls, text):
        """Return numbered anchors in `text` that are not a sequence starting at one.

        An orphan anchor refers to an item of a list that is no longer in the file,
        which is what a batch text edit leaves behind.
        """
        found = sorted({int(m.group(1)) for m in cls.NUMBERED.finditer(text)})
        if found and found != list(range(1, len(found) + 1)):
            return found
        return []

    @classmethod
    def anchor_lines(cls, path):
        return [no for no, line in cls._prose_lines(path) if cls.ANCHOR.search(line)]

    @classmethod
    def empty_lines(cls, path):
        return [no for no, line in cls._prose_lines(path) if cls.EMPTY.match(line.strip())]

    @classmethod
    def stray_anchor_lines(cls, path):
        return [no for no, line in cls._prose_lines(path) if cls.STRAY_ANCHOR.match(line.strip())]

    @classmethod
    def open_paren_lines(cls, path):
        return [no for no, line in cls._prose_lines(path) if cls.OPEN_PAREN.search(line)]

    @classmethod
    def unbalanced_units(cls, path):
        return [
            start
            for start, text in cls._prose_units(path)
            if cls.paren_counts(text)[0] != cls.paren_counts(text)[1]
        ]

    @classmethod
    def truncated_units(cls, path):
        return [
            (start, text[-60:])
            for start, text in cls._prose_units(path)
            if cls.DANGLING.search(text)
        ]

    @classmethod
    def stray_anchor_units(cls, path):
        return [
            (start, cls.stray_numbered_anchors(text))
            for start, text in cls._prose_units(path)
            if cls.stray_numbered_anchors(text)
        ]

    @classmethod
    def sentence_tail_anchor_units(cls, path):
        return [
            (start, text[max(0, m.start() - 40) : m.end()])
            for start, text in cls._prose_units(path)
            for m in [cls.TAIL_ANCHOR.search(text)]
            if m
        ]

    # ---------------- the scans themselves ----------------

    # One case per rule rather than one per rule per file: the assertion names
    # every offending file and the line inside it, so a failure is no harder to
    # locate than under a per-file parametrisation, and the collected count stays
    # a statement about distinct behaviours checked.

    @classmethod
    def _scan(cls, predicate, message):
        """Run `predicate` over every scanned file and report each hit in place."""
        bad = []
        for path in cls.TARGETS:
            found = predicate(path)
            if found:
                bad.append("%s: %s at %s" % (path.relative_to(REPO), message, found))
        assert not bad, "\n  ".join(bad)

    def test_no_space_before_sentence_period(self):
        self._scan(self.anchor_lines, "stray space before a sentence period")

    def test_no_orphan_period_line(self):
        self._scan(self.empty_lines, "prose line that is only a period")

    def test_no_orphan_numbered_anchor_line(self):
        self._scan(self.stray_anchor_lines, "prose line that is only a numbered anchor")

    def test_no_space_after_an_opening_parenthesis(self):
        self._scan(self.open_paren_lines, "a word lost right after an opening parenthesis")

    def test_docstrings_and_comment_blocks_have_balanced_parentheses(self):
        self._scan(self.unbalanced_units, "unbalanced parentheses in a docstring or comment")

    def test_docstrings_and_comment_blocks_are_not_truncated(self):
        self._scan(self.truncated_units, "sentence cut off mid-phrase")

    def test_numbered_anchors_point_at_a_list(self):
        self._scan(self.stray_anchor_units, "numbered anchor with no matching list")

    def test_no_sentence_ends_on_a_numbered_anchor(self):
        self._scan(
            self.sentence_tail_anchor_units,
            "sentence that ends on a numbered anchor and points at nothing",
        )

    # ---------------- proof that each rule bites ----------------
    #
    # Samples are assembled from fragments, because this file is one of the files
    # the scan covers and a literal damaged sentence would trip its own guard.

    @staticmethod
    def _write(tmp_path, name, text):
        path = tmp_path / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_guard_sees_a_truncated_line_at_column_zero(self, tmp_path):
        """The old indentation floor hid module docstring lines: the class defect."""
        broken = self._write(
            tmp_path,
            "broken.py",
            '"""Summary.\n\nSilent fallback was the largest source of wrong XML in the\n"""\n'
            "x = 1\n",
        )
        fixed = self._write(
            tmp_path,
            "fixed.py",
            '"""Summary.\n\nSilent fallback was the largest source of wrong XML.\n"""\n' "x = 1\n",
        )
        assert [start for start, _ in self.truncated_units(broken)] == [1]
        assert not self.truncated_units(fixed)

    def test_guard_sees_an_unbalanced_parenthesis_in_a_docstring(self, tmp_path):
        """The old block scan covered comments only: the second class defect."""
        broken = self._write(
            tmp_path,
            "broken.py",
            'def f():\n    """Digests differ for two datasets'
            + " ("
            + ' /.\n\n    Returns:\n        Nothing.\n    """\n'
            "    return None\n",
        )
        fixed = self._write(
            tmp_path,
            "fixed.py",
            'def f():\n    """Digests differ.\n\n    Returns:\n        Nothing.\n    """\n'
            "    return None\n",
        )
        assert self.unbalanced_units(broken) == [2]
        assert not self.unbalanced_units(fixed)

    def test_guard_sees_orphan_punctuation(self, tmp_path):
        """A prose line that is only a period, or only a numbered anchor, was invisible."""
        broken = self._write(
            tmp_path,
            "broken.py",
            '"""Summary.\n.\n"""\n' + "# (3).\n" + "x = 1\n",
        )
        assert self.empty_lines(broken) == [2]
        assert self.stray_anchor_lines(broken) == [4]

    def test_guard_sees_a_stray_numbered_anchor_inside_prose(self, tmp_path):
        """An anchor with no list to point at is reported, a real list is not."""
        broken = self._write(
            tmp_path,
            "broken.py",
            "# declaration order must not decide whether the config is\n# valid (0).\nx = 1\n",
        )
        listed = self._write(
            tmp_path,
            "listed.py",
            "# Steps (1) and (2) both hold.\nx = 1\n",
        )
        assert self.stray_anchor_units(broken) == [(1, [0])]
        assert not self.stray_anchor_units(listed)

    def test_guard_sees_a_sentence_ending_on_its_own_anchor(self, tmp_path):
        """A merged list leaves one anchor behind, and one anchor is a valid sequence.

        The sequence rule above cannot report that shape, so the sentence-ending
        rule exists for it: a lone anchor closing a sentence refers to an item
        that is no longer in the file.
        """
        broken = self._write(
            tmp_path,
            "broken.py",
            '"""Read a number.\n\nThe guards were made unreachable' + " (1)" + '.\n"""\n',
        )
        listed = self._write(
            tmp_path,
            "listed.py",
            '"""Read a number.\n\nSteps (1) and (2) both apply to the guards.\n"""\n',
        )
        assert [start for start, _ in self.sentence_tail_anchor_units(broken)] == [1]
        assert not self.sentence_tail_anchor_units(listed)

    def test_guard_sees_a_word_lost_after_an_opening_parenthesis(self, tmp_path):
        """The alias comment lost its first word inside a parenthesis; no rule saw it."""
        broken = self._write(
            tmp_path,
            "broken.py",
            "# Public alias: the CLI and the library agree on what validated means"
            + " ("
            + " covered only the CLI before).\nx = 1\n",
        )
        fixed = self._write(
            tmp_path,
            "fixed.py",
            "# Public alias: the CLI and the library agree on what validated means"
            + " (the helper covered only the CLI before).\nx = 1\n",
        )
        assert self.open_paren_lines(broken) == [1]
        assert not self.open_paren_lines(fixed)

    def test_guard_sees_a_period_glued_to_a_space(self, tmp_path):
        """The one rule the old class had, re-proven so it cannot rot back."""
        broken = self._write(
            tmp_path,
            "broken.py",
            "def f():\n    " + '"""Every named choice follows the same rule' + " ." + '"""\n'
            "    return None\n",
        )
        assert self.anchor_lines(broken) == [2]

    def test_guard_scans_a_docstring_not_only_comments(self, tmp_path):
        """The first version never parsed docstrings at all."""
        broken = self._write(
            tmp_path,
            "broken.py",
            "class C:\n    "
            + '"""A Gamma beta must mean the same thing everywhere'
            + " ."
            + '"""\n'
            "    pass\n",
        )
        assert self.anchor_lines(broken) == [2]
        assert len(list(self._docstrings(broken))) == 1

    def test_the_target_list_is_not_silently_empty(self):
        """A guard parametrised over an empty list passes and proves nothing."""
        names = [str(p.relative_to(REPO)) for p in self.TARGETS]
        assert len(names) >= 30, names
        assert any(n.startswith("beast2py/") for n in names)
        assert any(n.startswith("tests/") for n in names)


class TestDistributionContents:
    """What the distributed release must contain, checked without building it.

    The wheel and sdist are assembled from MANIFEST.in plus the package-data table, so a
    licence file that is missing from those two declarations is missing from the
    archive users download, even though it is present in the working tree.
    """

    MANIFEST = (REPO / "MANIFEST.in").read_text(encoding="utf-8")
    PYPROJECT = (REPO / "pyproject.toml").read_text(encoding="utf-8")

    @classmethod
    def _included_paths(cls):
        """Every path an ``include`` line names, whether alone or beside others."""
        return {
            token
            for line in cls.MANIFEST.splitlines()
            if line.strip().startswith("include")
            for token in line.split()[1:]
        }

    def test_licence_and_notice_are_shipped(self):
        assert "beast2py/tools/NOTICE" in self.MANIFEST
        assert "LICENSE.BEAST2-LGPL-2.1.txt" in self.MANIFEST
        assert (PKG / "tools" / "NOTICE").exists()
        assert (PKG / "tools" / "LICENSE.BEAST2-LGPL-2.1.txt").exists()

    def test_example_fingerprint_pair_is_shipped(self):
        # the pairing is checkable "without running anything" only if both
        # files are in the archive the reader downloads.
        assert "examples/output_basic.xml" in self.MANIFEST
        assert "output_basic.fingerprint.json" in self.MANIFEST

    def test_both_validator_launchers_are_declared(self):
        """The archive must carry the launcher the test suite looks for by path.

        ``tests/test_beast2_validation.py`` resolves ``beast2_validate.sh`` at the
        repository root and skips all 21 BEAST2 integration tests when that one
        file is absent. The package copy alone is not enough: an unpacked sdist
        without the root wrapper therefore reports a fully green suite in which
        gate 3 never ran, and the skips are only visible under ``-rs``. The root
        file is a six-line wrapper that execs the package's own launcher, so both
        have to ship.
        """
        declared = self._included_paths()
        flat = " ".join(self.MANIFEST.split())
        assert "beast2_validate.sh" in declared, sorted(declared)
        assert "beast2py/beast2_validate.sh" in flat
        assert (REPO / "beast2_validate.sh").is_file()
        # and the package-data table must keep naming the inner one, which is
        # what puts the launcher into the wheel
        assert "beast2_validate.sh" in self.PYPROJECT

    def test_the_test_harness_is_declared_for_the_sdist(self):
        """A source install has to be able to *collect* the suite it advertises.

        setuptools matches the ``test_*`` modules but nothing else in ``tests/``,
        so an unpacked sdist shipped the tests without their fixtures and without
        the package marker: every module imports them relatively, so pytest aborted
        the whole run with "attempted relative import with no known parent
        package", which is a louder but no better outcome than the silent skip.
        """
        declared = self._included_paths()
        for harness in ("tests/conftest.py", "tests/__init__.py"):
            assert harness in declared, sorted(declared)

    def test_declared_license_is_the_same_expression_everywhere(self):
        """One SPDX expression, asserted at every place that declares a licence.

        The archived release declared its licence three different ways: the
        package metadata and the citation file said ``MIT AND LGPL-2.1-only``
        while ``.zenodo.json`` said bare ``MIT`` and ``beast2py.__license__`` said
        bare ``MIT``. The bare declarations understated the vendored LGPL
        component, and the archive record is the most authoritative place a reader
        looks. An earlier version of this check asserted only that the letters
        "LGPL" appeared somewhere in ``pyproject.toml``, which the SPDX string
        already satisfies -- an assertion that cannot fail is not an assertion.

        PEP 639 is why there is no ``License ::`` trove classifier: the SPDX
        expression replaces it, and the packaging job in ``ci.yml`` asserts the
        classifiers stay absent. Naming the licence twice, in two formats, is
        exactly how a three-way disagreement becomes possible.
        """
        combined = "MIT AND LGPL-2.1-only"
        sources = {
            "pyproject.toml": self.PYPROJECT,
            "CITATION.cff": (REPO / "CITATION.cff").read_text(encoding="utf-8"),
            ".zenodo.json": (REPO / ".zenodo.json").read_text(encoding="utf-8"),
            "README.md": (REPO / "README.md").read_text(encoding="utf-8"),
            "README_zh.md": (REPO / "README_zh.md").read_text(encoding="utf-8"),
            "beast2py/__init__.py": (PKG / "__init__.py").read_text(encoding="utf-8"),
            "LICENSE": (REPO / "LICENSE").read_text(encoding="utf-8"),
            "beast2py/tools/NOTICE": (PKG / "tools" / "NOTICE").read_text(encoding="utf-8"),
        }
        for label, text in sources.items():
            assert combined in text, f"{label} does not declare {combined}"

        # the bare declaration is the defect this closes: no place may state only
        # MIT, and "MIT" may not appear as a licence value without the LGPL half.
        assert '"license": {"id": "MIT"}' not in sources[".zenodo.json"]
        assert '__license__ = "MIT"' not in sources["beast2py/__init__.py"]
        # read the parsed value, not the raw text: pyproject.toml explains in a comment why
        # the License :: trove classifiers were removed, and a substring test would trip on
        # that explanation instead of on an actual classifier
        import tomllib

        classifiers = tomllib.loads(self.PYPROJECT).get("project", {}).get("classifiers", [])
        assert not any(
            c.startswith("License ::") for c in classifiers
        ), f"a License :: classifier must not sit beside the PEP 639 SPDX expression: {classifiers}"

    def test_the_bundled_classes_keep_their_own_licence_text(self):
        """LGPL-2.1 is only honoured if the text and the notice travel with it."""
        lgpl = PKG / "tools" / "LICENSE.BEAST2-LGPL-2.1.txt"
        assert lgpl.is_file(), "the vendored BEAST2 classes lost their licence text"
        head = lgpl.read_text(encoding="utf-8", errors="replace")
        assert "GNU LESSER GENERAL PUBLIC LICENSE" in head
        assert "Version 2.1" in head
        notice = (PKG / "tools" / "NOTICE").read_text(encoding="utf-8")
        assert "beast/pkgmgmt" in notice
        assert "LGPL" in notice and "2.1" in notice
        # and both have to be declared, or they are in the tree but not in the
        # archive a reader downloads
        assert "beast2py/tools/NOTICE" in self.MANIFEST
        assert "LICENSE.BEAST2-LGPL-2.1.txt" in self.MANIFEST


class TestRepositoryHygiene:
    """What the repository publishes is its public source.

    Two failure modes are cheap to prevent and expensive to retract once
    published: an absolute path from the author's machine, and a
    credential committed by accident.  The scan covers exactly the files git
    would ship - tracked plus untracked-but-not-ignored - so agent state,
    build residue and caches are out of scope by construction, not by a
    hand-maintained exclusion list.
    """

    SUFFIXES = {
        ".py",
        ".md",
        ".toml",
        ".txt",
        ".yaml",
        ".yml",
        ".sh",
        ".java",
        ".json",
        ".in",
        ".cfg",
        ".xml",
    }
    SECRET = re.compile(
        r"(?i)\b(api[_-]?key|secret[_-]?key|access[_-]?token|auth[_-]?token|"
        r"client[_-]?secret|private[_-]?key|password)\b\s*[:=]\s*\S"
        r"|\bgh[pousr]_[A-Za-z0-9]{18,}\b|\bAKIA[0-9A-Z]{16}\b"
        r"|-----BEGIN [A-Z ]*PRIVATE KEY-----|\bxox[baprs]-[A-Za-z0-9-]{10,}"
    )
    HOME = re.compile(
        r"(?<!\w)(?:/Users/[\w.\-]+/|/home/[\w.\-]+/" r"|[A-Za-z]:\\\\Users\\\\[\w.\-]+\\\\)"
    )

    @staticmethod
    def _shipped_text_files():
        import subprocess

        try:
            listed = subprocess.run(
                ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                cwd=str(REPO),
                capture_output=True,
                text=True,
                check=True,
            )
        except Exception:  # no git: skip cleanly
            pytest.skip("git is not available to enumerate shipped files")
        names = [
            n
            for n in listed.stdout.splitlines()
            if pathlib.Path(n).suffix.lower() in TestRepositoryHygiene.SUFFIXES
            and not pathlib.Path(n).name.startswith("LICENSE.BEAST2")
        ]
        return [REPO / n for n in names if (REPO / n).is_file()]

    def test_the_scan_is_not_silently_empty(self):
        files = self._shipped_text_files()
        assert len(files) > 30, "expected to scan the package and docs, found %d" % len(files)

    def test_no_absolute_paths_from_a_machines(self):
        hits = []
        for path in self._shipped_text_files():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for m in self.HOME.finditer(text):
                ctx = text[max(0, m.start() - 40) : m.end() + 30].replace("\n", " ")
                hits.append("%s: %s" % (path.relative_to(REPO), ctx))
        assert not hits, "absolute home paths in shipped files:\n  %s" % "\n  ".join(hits[:8])

    def test_no_credential_shaped_strings(self):
        hits = []
        for path in self._shipped_text_files():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for m in self.SECRET.finditer(text):
                hits.append("%s: %s" % (path.relative_to(REPO), m.group(0)[:60]))
        assert not hits, "credential-shaped strings:\n  %s" % "\n  ".join(hits[:8])

    def test_the_patterns_reject_what_they_claim_to(self):
        # A guard that cannot fail is worse than no guard, so prove it bites.
        # The samples are assembled from fragments so that this file, which the
        # scan also covers, does not report itself.
        home = "/".join(["", "Users", "someone", "beast2py", "x.py"])
        lin = "/".join(["", "home", "someone", "beast2py", "y"])
        cred = "api" + "_" + "key = " + "ghp_" + "a" * 20
        pw = "pass" + "word: hunter2"
        assert self.HOME.search(home), home
        assert self.HOME.search(lin), lin
        assert self.SECRET.search(cred)
        assert self.SECRET.search(pw)
        assert not self.SECRET.search("tokens = value.split()")
        assert not self.SECRET.search("a change-detection token, not a key")

    def test_no_junk_is_tracked(self):
        import subprocess

        tracked = subprocess.run(
            ["git", "ls-files"], cwd=str(REPO), capture_output=True, text=True
        ).stdout.splitlines()
        junk = [
            n
            for n in tracked
            if pathlib.Path(n).name in {".DS_Store"} or n.endswith((".pyc", ".mimosa"))
        ]
        assert not junk, "junk files are tracked: %s" % junk[:5]


class TestTutorialWalkthrough:
    """The tutorials are numbered steps, so they have to survive being run as one script."""

    TUTORIALS = ["docs/tutorial_en.md", "docs/tutorial_zh.md"]
    FENCE = re.compile(r"```[a-zA-Z_+-]*\n(.*?)```", re.S)
    # `--xml` names an input to validate; `--output-dir` is a directory that the
    # three documented --type variants legitimately share.
    WRITERS = ("--output", "-o", "--report")

    def _commands(self, rel):
        text = (REPO / rel).read_text(encoding="utf-8")
        cmds = []
        for block in self.FENCE.finditer(text):
            for line in block.group(1).replace("\\\n", " ").split("\n"):
                line = re.sub(r"\s+", " ", line.strip())
                if not line.startswith("beast2py "):
                    continue
                if re.search(r"\[[^\]]*\]", line[len("beast2py") :]):
                    continue  # a synopsis line, not something to run
                if "--help" in line or "--version" in line:
                    continue
                cmds.append(line)
        return cmds

    def _writes(self, cmd):
        toks = cmd.split()
        return [toks[i + 1] for i, t in enumerate(toks) if t in self.WRITERS and i + 1 < len(toks)]

    @pytest.mark.parametrize("rel", TUTORIALS)
    def test_no_two_steps_write_the_same_file(self, rel):
        # generate/quick refuse to replace an output without --force, so two
        # steps sharing a name would leave the reader at a refusal.
        written = {}
        for cmd in self._commands(rel):
            if "--force" in cmd.split():
                continue
            for path in self._writes(cmd):
                written.setdefault(path, []).append(cmd)
        clash = {p: c for p, c in written.items() if len(c) > 1}
        assert not clash, "%s writes %s more than once: %s" % (
            rel,
            sorted(clash),
            clash[sorted(clash)[0]][0],
        )

    @pytest.mark.parametrize("rel", TUTORIALS)
    def test_the_walkthrough_names_commands_that_exist(self, rel):
        from beast2py.main import build_parser

        parser = build_parser()
        choices = next(a.choices for a in parser._actions if a.dest == "command")
        for cmd in self._commands(rel):
            name = cmd.split()[1]
            assert name in choices, "%s documents unknown command: %s" % (rel, cmd)


class TestErrorDocumentation:
    """Troubleshooting section 22 promises the messages the code actually emits.

    A table that drifts from the code is worse than no table: the reader trusts it to
    diagnose a failure, searches for the documented wording, and finds no match.
    """

    MESSAGES = (
        "Invalid YAML in",
        "Calibration file not found",
        "must be a list of calibration points",
    )
    MANUALS = [
        ("docs/user_manual_en.md", "## 22. Troubleshooting"),
        ("docs/user_manual_zh.md", "## 22. \u6545\u969c\u6392\u9664"),
    ]

    def test_the_code_still_emits_these_messages(self, tmp_path):
        from beast2py.config import ConfigError, ConfigParser
        from beast2py.main import _parse_calibration_yaml

        bad = tmp_path / "bad.yaml"
        bad.write_text("a: [unclosed\n", encoding="utf-8")
        with pytest.raises(ConfigError) as cfg:
            ConfigParser().parse(str(bad))
        assert "Invalid YAML in" in str(cfg.value)
        assert str(bad) in str(cfg.value), "the message must name the file"

        with pytest.raises(ConfigError) as miss:
            _parse_calibration_yaml(str(tmp_path / "nope.yaml"))
        assert "Calibration file not found" in str(miss.value)

        mapping = tmp_path / "map.yaml"
        mapping.write_text("name: x\n", encoding="utf-8")
        with pytest.raises(ConfigError) as shape:
            _parse_calibration_yaml(str(mapping))
        assert "must be a list" in str(shape.value)

    @pytest.mark.parametrize("rel,heading", MANUALS)
    def test_both_manuals_document_them(self, rel, heading):
        text = (REPO / rel).read_text(encoding="utf-8")
        assert heading in text, "%s has no troubleshooting section" % rel
        section = text[text.index(heading) :]
        missing = [m for m in self.MESSAGES if m not in section]
        assert not missing, "%s does not document: %s" % (rel, missing)


class TestPublicSurface:
    """What a reader may import, and with which defaults, is itself a claim.

    The article advertises the CLI "and equivalent Python functions". Two places made that
    claim weaker than it reads: the package root did not re-export the data classes the API
    returns, so `from beast2py import Alignment` raised ImportError, and `read_sequence`
    defaulted to declared nucleotide data where `beast2py quick` defaults to inferring the
    type from the residues -- an amino-acid alignment worked from the CLI and silently did not
    work from the function of the same name. `api.py` also had no `__all__`, so nothing in the
    package said which of its module-level names were public at all.
    """

    ROOT_REEXPORTS = ["Alignment", "Sequence", "HyperpriorConfig", "NON_SCIENCE_METADATA"]

    def test_package_root_reexports_the_data_classes_the_api_returns(self):
        import beast2py.api as api

        aln = api.read_sequence(str(EXAMPLES / "primates.fasta"))
        assert isinstance(aln, beast2py.Alignment), "the returned class is not importable by name"
        for name in self.ROOT_REEXPORTS:
            assert hasattr(beast2py, name), "from beast2py import %s fails" % name
            assert name in beast2py.__all__, "%s is importable but not declared public" % name

    def test_both_public_surfaces_are_declared_and_honest(self):
        import beast2py.api as api

        for mod in (beast2py, api):
            promised = [n for n in mod.__all__ if not hasattr(mod, n)]
            assert not promised, "%s.__all__ promises %s that the module lacks" % (
                mod.__name__,
                promised,
            )

    def test_api_declares_every_public_name_it_defines(self):
        """A later addition to api.py must be declared public or this test fails."""
        import ast
        import beast2py.api as api

        tree = ast.parse((PKG / "api.py").read_text(encoding="utf-8"))
        defined = sorted(
            n.name
            for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.ClassDef)) and not n.name.startswith("_")
        )
        undeclared = [n for n in defined if n not in api.__all__]
        assert not undeclared, "api.py defines %s but keeps them off __all__" % undeclared

    def test_read_sequence_defaults_to_the_same_inference_as_the_cli(self):
        import inspect
        import beast2py.api as api
        from beast2py.main import build_parser
        from beast2py.models import DataType

        cli_default = (
            build_parser().parse_args(["quick", "-a", "in.fasta", "-o", "out.xml"]).data_type
        )
        assert cli_default == "auto", "the CLI stopped defaulting to inference"
        for target in (api.read_sequence, api.Beast2Py.read_sequence):
            default = inspect.signature(target).parameters["data_type"].default
            assert default == cli_default, "%s defaults to %r, the CLI defaults to %r" % (
                target.__qualname__,
                default,
                cli_default,
            )

        inferred = api.read_sequence(str(EXAMPLES / "aminoacid.fasta"))
        assert inferred.data_type is DataType.AMINOACID, (
            "reading the shipped amino-acid example without arguing about its type is the "
            "behaviour this parity claim is about; got %r" % (inferred.data_type,)
        )
        nucleotides = api.read_sequence(str(EXAMPLES / "primates.fasta"))
        assert nucleotides.data_type is DataType.NUCLEOTIDE

    def test_quick_mode_is_the_same_function_from_either_entry_point(self, tmp_path):
        """The API's quick mode used to read every alignment as nucleotides.

        ``beast2py quick`` infers the data type unless told, which is what its own help text
        says makes a protein alignment work at all. The Python function of the same name had no
        such knob and hardcoded a declared nucleotide type, so the article's "equivalent Python
        functions" was false for the shipped amino-acid example: the two entry points emitted
        different XML from the same input.
        """
        from beast2py.api import Beast2Py
        from beast2py.main import main

        aa = str(EXAMPLES / "aminoacid.fasta")
        cli_out = tmp_path / "cli.xml"
        api_out = tmp_path / "api.xml"
        assert (
            main(
                [
                    "quick",
                    "-a",
                    aa,
                    "-o",
                    str(cli_out),
                    "--subst-model",
                    "jtt",
                    "--clock-model",
                    "strict",
                ]
            )
            == 0
        ), "quick mode on the shipped protein example must work from the CLI"
        Beast2Py().quick_generate(
            alignment=aa,
            output=str(api_out),
            subst_model="jtt",
            clock_model="strict",
        )
        cli_xml = cli_out.read_text(encoding="utf-8")
        api_xml = api_out.read_text(encoding="utf-8")
        assert cli_xml == api_xml, (
            "quick mode must not depend on which entry point was used; the CLI emits %d "
            "bytes, the API %d" % (len(cli_xml), len(api_xml))
        )
        assert "Dayhoff" in cli_xml or "JTT" in cli_xml, (
            "the protein model must be present in either output, so this comparison cannot "
            "pass by both sides being equally wrong"
        )
