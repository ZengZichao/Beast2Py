"""Release-integrity guards derived from the submission review.

Each test locks a defect that the review found, so it cannot silently return:

* the birthday-collision bound quoted to users (ISSUE-001),
* comment text that lost a character to an earlier edit pass (ISSUE-023),
* the sidecar/XML pair that lets a reader verify an analysis without running
  anything (ISSUE-031),
* order-stability of the fingerprint serialiser (ISSUE-030),
* what the built wheel and sdist must contain, checked from the declarations
  (ISSUE-004, 015, 031),
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
        space = 2 ** bits
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
        than the committed example the paper says was validated.
        """
        import contextlib
        import io
        import tempfile

        from beast2py.api import generate_xml
        from beast2py.config import ConfigParser

        out = pathlib.Path(tempfile.mkdtemp()) / "output_basic.xml"
        with contextlib.redirect_stdout(io.StringIO()):
            generate_xml(str(EXAMPLES / "config_basic.yaml"), str(out), force=True)
        assert out.read_bytes() == (EXAMPLES / "output_basic.xml").read_bytes(), (
            "the committed exemplar XML is behind the writer")

        config = ConfigParser().parse(EXAMPLES / "config_basic.yaml")
        side = json.loads((EXAMPLES / "output_basic.fingerprint.json").read_text())
        fresh = FingerprintGenerator.generate_fingerprint_dict(
            config, xml_content=out.read_text(encoding="utf-8"))
        volatile = {"generation_time"}
        assert {k: v for k, v in side.items() if k not in volatile} == \
            {k: v for k, v in fresh.items() if k not in volatile}, (
            "the committed sidecar is behind what the shipped code writes")

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

        one = json.dumps({"b": 1, "a": 2}, default=FingerprintGenerator._canonical,
                         sort_keys=True)
        two = json.dumps({"a": 2, "b": 1}, default=FingerprintGenerator._canonical,
                         sort_keys=True)
        assert one == two


class TestCommentIntegrity:
    """Earlier edit passes stripped characters from comments and docstrings.

    The failures were mechanical and therefore are mechanically detectable, so
    they are pinned here rather than left to review.
    """

    # A period glued to nothing but preceded by a space is the signature of an
    # edit that deleted the word after it.  The alternatives are what a repaired
    # sentence can legitimately end with, and deliberately exclude "from . import".
    ANCHOR = re.compile(r"(\w|[^\x00-\x7f]) +\.(?:$| [A-Z`#]| [-\u2014])")
    EMPTY = re.compile(r"#\s*\.\s*$")

    def _prose_lines(self, path):
        for no, line in enumerate(path.read_text(encoding="utf-8").split("\n"), 1):
            stripped = line.strip()
            if stripped.startswith("#") or re.match(r"^\s{4,}\S", line):
                yield no, line

    @pytest.mark.parametrize(
        "path", sorted(PKG.rglob("*.py")), ids=lambda p: p.name
    )
    def test_no_space_before_sentence_period(self, path):
        bad = [n for n, line in self._prose_lines(path) if self.ANCHOR.search(line)]
        assert not bad, f"{path.name}: stray space before '.' at line(s) {bad}"

    @pytest.mark.parametrize(
        "path", sorted(PKG.rglob("*.py")), ids=lambda p: p.name
    )
    def test_no_orphan_period_comment(self, path):
        bad = [
            n
            for n, line in self._prose_lines(path)
            if self.EMPTY.match(line.strip())
        ]
        assert not bad, f"{path.name}: comment that is only a '.' at line(s) {bad}"

    @pytest.mark.parametrize(
        "path", sorted(PKG.rglob("*.py")), ids=lambda p: p.name
    )
    def test_comment_blocks_have_balanced_parentheses(self, path):
        lines = path.read_text(encoding="utf-8").split("\n")
        block, start, problems = [], 0, []
        for no, line in enumerate(lines, 1):
            if line.strip().startswith("#"):
                if not block:
                    start = no
                block.append(line.strip().lstrip("# "))
                continue
            if block:
                text = " ".join(block)
                # interval notation such as "(0, 1]" legitimately unbalances.
                cleaned = re.sub(r"\(\s*[\d.,\s]+?[\]\)]", "", text)
                if cleaned.count("(") != cleaned.count(")"):
                    problems.append((start, cleaned.count("("), cleaned.count(")")))
                block = []
        assert not problems, f"{path.name}: unbalanced parentheses at {problems}"


class TestDistributionContents:
    """What the deposited release must contain, checked without building it.

    The wheel and sdist are assembled from MANIFEST.in plus the package-data table, so a
    licence file or authoring gate that is missing from those two declarations is missing
    from the archive the paper cites, even though it is present in the working tree.
    """

    MANIFEST = (REPO / "MANIFEST.in").read_text(encoding="utf-8")
    PYPROJECT = (REPO / "pyproject.toml").read_text(encoding="utf-8")

    def test_licence_and_notice_are_shipped(self):
        assert "beast2py/tools/NOTICE" in self.MANIFEST
        assert "LICENSE.BEAST2-LGPL-2.1.txt" in self.MANIFEST
        assert (PKG / "tools" / "NOTICE").exists()
        assert (PKG / "tools" / "LICENSE.BEAST2-LGPL-2.1.txt").exists()

    def test_authoring_gates_are_shipped(self):
        # ISSUE-015: the draw.io export gate has to travel with the deposited version.
        assert "recursive-include scripts *.py" in self.MANIFEST
        for script in ("check_figure_export.py", "check_submission_placeholders.py"):
            assert (REPO / "scripts" / script).exists(), script

    def test_example_fingerprint_pair_is_shipped(self):
        # ISSUE-031: the pairing is checkable "without running anything" only if both
        # files are in the archive the reader downloads.
        assert "examples/output_basic.xml" in self.MANIFEST
        assert "output_basic.fingerprint.json" in self.MANIFEST

    def test_declared_license_matches_the_bundled_classes(self):
        assert "MIT AND LGPL-2.1-only" in self.PYPROJECT
        assert "LGPL" in self.PYPROJECT  # the classifier, not only the SPDX string


class TestRepositoryHygiene:
    """What the repository publishes is the archived source of the paper.

    Two failure modes are cheap to prevent and expensive to retract once a
    deposit DOI is minted: an absolute path from the author's machine, and a
    credential committed by accident.  The scan covers exactly the files git
    would ship - tracked plus untracked-but-not-ignored - so agent state,
    build residue and caches are out of scope by construction, not by a
    hand-maintained exclusion list.
    """

    SUFFIXES = {".py", ".md", ".toml", ".txt", ".yaml", ".yml", ".sh", ".java",
                ".json", ".in", ".cfg", ".xml"}
    SECRET = re.compile(
        r"(?i)\b(api[_-]?key|secret[_-]?key|access[_-]?token|auth[_-]?token|"
        r"client[_-]?secret|private[_-]?key|password)\b\s*[:=]\s*\S"
        r"|\bgh[pousr]_[A-Za-z0-9]{18,}\b|\bAKIA[0-9A-Z]{16}\b"
        r"|-----BEGIN [A-Z ]*PRIVATE KEY-----|\bxox[baprs]-[A-Za-z0-9-]{10,}")
    HOME = re.compile(r"(?<!\w)(?:/Users/[\w.\-]+/|/home/[\w.\-]+/"
                      r"|[A-Za-z]:\\\\Users\\\\[\w.\-]+\\\\)")

    @staticmethod
    def _shipped_text_files():
        import subprocess
        try:
            listed = subprocess.run(
                ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
                cwd=str(REPO), capture_output=True, text=True, check=True)
        except Exception:                                  # no git: skip cleanly
            pytest.skip("git is not available to enumerate shipped files")
        names = [n for n in listed.stdout.splitlines()
                 if pathlib.Path(n).suffix.lower() in TestRepositoryHygiene.SUFFIXES
                 and not pathlib.Path(n).name.startswith("LICENSE.BEAST2")]
        return [REPO / n for n in names if (REPO / n).is_file()]

    def test_the_scan_is_not_silently_empty(self):
        files = self._shipped_text_files()
        assert len(files) > 30, "expected to scan the package and docs, found %d" % len(files)

    def test_no_absolute_paths_from_a_machines(self):
        hits = []
        for path in self._shipped_text_files():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for m in self.HOME.finditer(text):
                ctx = text[max(0, m.start() - 40):m.end() + 30].replace("\n", " ")
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
        tracked = subprocess.run(["git", "ls-files"], cwd=str(REPO),
                                 capture_output=True, text=True).stdout.splitlines()
        junk = [n for n in tracked
                if pathlib.Path(n).name in {".DS_Store"} or n.endswith((".pyc", ".mimosa"))]
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
                if re.search(r"\[[^\]]*\]", line[len("beast2py"):]):
                    continue  # a synopsis line, not something to run
                if "--help" in line or "--version" in line:
                    continue
                cmds.append(line)
        return cmds

    def _writes(self, cmd):
        toks = cmd.split()
        return [toks[i + 1] for i, t in enumerate(toks)
                if t in self.WRITERS and i + 1 < len(toks)]

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
            rel, sorted(clash), clash[sorted(clash)[0]][0])

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

    MESSAGES = ("Invalid YAML in",
                "Calibration file not found",
                "must be a list of calibration points")
    MANUALS = [("docs/user_manual_en.md", "## 22. Troubleshooting"),
               ("docs/user_manual_zh.md", "## 22. \u6545\u969c\u6392\u9664")]

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
        section = text[text.index(heading):]
        missing = [m for m in self.MESSAGES if m not in section]
        assert not missing, "%s does not document: %s" % (rel, missing)
