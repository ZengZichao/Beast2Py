"""XML validation module for Beast2Py.

Validates generated XML for structural correctness, id/idref consistency,
and BEAST2-specific requirements.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from defusedxml import ElementTree as SafeElementTree
from defusedxml.common import DefusedXmlException
from pathlib import Path
from typing import List, Optional, Tuple


class ValidationResult:
    """Container for validation results."""

    def __init__(self):
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def add_error(self, msg: str) -> None:
        self.errors.append(msg)

    def add_warning(self, msg: str) -> None:
        self.warnings.append(msg)

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0

    def __str__(self) -> str:
        lines = []
        if self.errors:
            lines.append(f"ERRORS ({len(self.errors)}):")
            for e in self.errors:
                lines.append(f"  ✗ {e}")
        if self.warnings:
            lines.append(f"WARNINGS ({len(self.warnings)}):")
            for w in self.warnings:
                lines.append(f"  ⚠ {w}")
        if not self.errors and not self.warnings:
            lines.append("✓ XML is valid. No errors or warnings.")
        return "\n".join(lines)


class XMLValidator:
    """Validate BEAST2 XML files."""

    @staticmethod
    def validate(xml_string: str, config=None) -> ValidationResult:
        """Validate an XML string for BEAST2 compatibility.

        Args:
            xml_string: The XML string to validate.
            config: Optional BEASTConfig, accepted for call-site symmetry; the
                semantic checks in :meth:`_check_semantic` run either way.

        Returns:
            ValidationResult with errors and warnings.
        """
        result = ValidationResult()

        # 1. XML format validation (defused parser: entity expansion rejected)
        try:
            root = SafeElementTree.fromstring(xml_string)
        except ET.ParseError as e:
            result.add_error(f"XML parse error: {e}")
            return result
        except DefusedXmlException as e:
            result.add_error(f"XML entity expansion is not allowed: {e}")
            return result

        # 2. Root element check
        if root.tag != "beast":
            result.add_error(f"Root element must be <beast>, found <{root.tag}>")
            return result

        # 3. Check for required elements
        has_data = False
        has_run = False
        for child in root:
            if child.tag == "data":
                has_data = True
            elif child.tag == "run":
                has_run = True
            elif child.tag == "distribution":
                # Posterior might be defined inline
                pass

        if not has_data:
            result.add_error("No <data> element found")
        if not has_run:
            result.add_error("No <run> element found")

        # 4. Collect all IDs and idrefs
        all_ids = set()
        all_idrefs = set()

        XMLValidator._collect_ids(root, all_ids, all_idrefs)

        # 4a. Collect taxon names from sequence elements (implicit IDs)
        for data_elem in root.findall("data"):
            for seq_elem in data_elem.findall("sequence"):
                taxon = seq_elem.get("taxon")
                if taxon:
                    all_ids.add(taxon)

        # 5. Check idref consistency
        for idref in all_idrefs:
            if idref not in all_ids:
                result.add_error(f"idref '{idref}' does not reference any defined id")

        # 6. Check for duplicate IDs
        id_counts: dict[str, int] = {}
        XMLValidator._count_ids(root, id_counts)
        for id_val, count in id_counts.items():
            if count > 1:
                result.add_error(f"Duplicate id: '{id_val}' appears {count} times")

        # 7. Check run element structure
        run_elem = root.find("run")
        if run_elem is not None:
            # Check for state
            state = run_elem.find("state")
            if state is None:
                result.add_error("No <state> element in <run>")

            # Check for distribution reference
            has_dist_ref = False
            for child in run_elem:
                if child.tag == "distribution" and "idref" in child.attrib:
                    has_dist_ref = True
                    break
            if not has_dist_ref:
                result.add_error("No posterior distribution reference in <run>")

            # Check for at least one logger
            loggers = run_elem.findall("logger")
            if len(loggers) == 0:
                result.add_error("No <logger> elements in <run>")

        # 8. Check data element has sequences
        for data_elem in root.findall("data"):
            seqs = data_elem.findall("sequence")
            if (
                len(seqs) == 0
                and data_elem.get("spec") != "beast.base.evolution.alignment.FilteredAlignment"
            ):
                result.add_warning(f"data '{data_elem.get('id')}' has no sequences")

        # 9. Check namespace attribute
        ns = root.get("namespace", "")
        if not ns:
            result.add_warning("No namespace attribute on <beast> element")
        elif "beast.base" not in ns:
            result.add_warning(
                f"Namespace does not contain 'beast.base' (may be v2.6.x format): {ns[:50]}..."
            )

        # 10. Semantic checks: coverage of priors, alignment and bound sanity
        XMLValidator._check_semantic(root, result)

        return result

    @staticmethod
    def _check_semantic(root: ET.Element, result: "ValidationResult") -> None:
        """Check the invariants that make an analysis *correct*, not just parseable.

        Every check here exists because a model was found that BEAST2
        happily parsed while computing something other than what the user
        specified.

        Args:
            root: The parsed ``<beast>`` element.
            result: Validation result to populate.
        """
        state = None
        for run_elem in root.findall("run"):
            state = run_elem.find("state")
            if state is not None:
                break
        if state is None:
            return

        # Estimated state nodes and the ids covered by a Prior
        estimated = []
        for elem in state:
            if elem.tag in ("parameter", "tree") and elem.get("estimate", "true") == "true":
                node_id = elem.get("id")
                if node_id:
                    estimated.append((node_id, elem))

        covered = set()
        for prior_elem in root.iter("distribution"):
            # Match the *simple* class name, case-insensitively. Only accepting
            # the fully qualified ".Prior" spelling made `beast2py validate`
            # report five invented "has no <Prior>" errors on its own example
            # once the file was re-spelled the way the BEAST2 GUI and the
            # official templates spell it.
            spec = prior_elem.get("spec", "")
            if spec.rsplit(".", 1)[-1].lower() != "prior":
                continue
            target = (prior_elem.get("x") or "").lstrip("@")
            if not target:
                # Nested form: <x name="input" idref="someNode"/>
                for child in prior_elem:
                    if child.get("name") in ("x", "input"):
                        target = child.get("idref") or (child.get("@ref") or "").lstrip("@")
                        if target:
                            break
            if target:
                covered.add(target)

        for node_id, elem in estimated:
            if elem.tag == "tree":
                continue  # the tree's density comes from the tree prior
            spec = elem.get("spec", "")
            if spec.endswith("IntegerParameter") or spec.endswith("BooleanParameter"):
                # Discrete parameters with a bounded support are uniform by
                # construction over their domain: rate category indices and
                # indicator bits need no explicit density.
                continue
            if node_id not in covered:
                result.add_error(
                    f"Estimated parameter '{node_id}' has no <Prior>: a flat prior on an "
                    f"unbounded support is improper, and the analysis cannot be "
                    f"reproduced from the XML because its prior is not written down"
                )

        for node_id, elem in estimated:
            lower = elem.get("lower")
            upper = elem.get("upper")
            if lower is not None and upper is not None:
                try:
                    if float(lower) >= float(upper):
                        result.add_error(
                            f"State node '{node_id}': lower ({lower}) >= upper ({upper})"
                        )
                except ValueError:
                    result.add_error(f"State node '{node_id}': non-numeric bounds")
            value = elem.get("value")
            if value is None:
                continue
            tokens = [tok for tok in value.replace(",", " ").split()]
            # dimension="auto" is legal BEAST2, and an unguarded int() here
            # crashed the whole validator instead of reporting anything.
            raw_dimension = elem.get("dimension", "1")
            try:
                dimension = int(raw_dimension)
            except ValueError:
                dimension = None
                if raw_dimension.strip().lower() not in ("auto", ""):
                    result.add_error(
                        f"State node '{node_id}': dimension must be a positive "
                        f"integer or 'auto', got '{raw_dimension}'"
                    )
            if dimension is None:
                continue
            try:
                numeric = [float(tok) for tok in tokens]
            except ValueError:
                numeric = []  # boolean / textual initial values
            if numeric and len(numeric) == 1 and dimension > 1:
                numeric = numeric * dimension
            if numeric and len(numeric) != dimension:
                result.add_error(
                    f"State node '{node_id}': {len(numeric)} initial values but "
                    f"dimension={dimension}"
                )
            for num in numeric:
                if lower is not None and num < float(lower):
                    result.add_error(
                        f"State node '{node_id}': initial value {num} is below lower ({lower})"
                    )
                    break
                if upper is not None and num > float(upper):
                    result.add_error(
                        f"State node '{node_id}': initial value {num} is above upper ({upper})"
                    )
                    break

        # Alignment sanity inside the produced XML ( backstop)
        for data_elem in root.findall("data"):
            seqs = data_elem.findall("sequence")
            if not seqs:
                continue
            lengths = {len((seq.text or "").strip().replace(" ", "")) for seq in seqs}
            if len(lengths) > 1:
                result.add_error(
                    f"<data id='{data_elem.get('id')}': sequences of unequal length "
                    f"{sorted(lengths)}"
                )
            taxa = [seq.get("taxon") for seq in seqs]
            dupes = sorted({t for t in taxa if taxa.count(t) > 1})
            if dupes:
                result.add_error(
                    f"<data id='{data_elem.get('id')}': duplicate taxa {', '.join(map(str, dupes))}"
                )

        # Resume capability
        for run_elem in root.findall("run"):
            store = run_elem.get("storeEvery")
            if store is None:
                continue
            try:
                store_value = int(store)
            except ValueError:
                result.add_error(
                    f"<run storeEvery='{store}'>: must be an integer"
                )
                continue
            if store_value == 0:
                result.add_error(
                    "<run storeEvery='0'>: 0 is neither a positive checkpoint "
                    "interval nor the -1 that disables checkpointing"
                )
            elif store_value < 0:
                result.add_warning(
                    f"<run storeEvery='{store}'>: state checkpointing is off, so this chain "
                    f"cannot be resumed after an interruption"
                )

    @staticmethod
    def validate_file(xml_file: str) -> ValidationResult:
        """Validate an XML file.

        Args:
            xml_file: Path to the XML file.

        Returns:
            ValidationResult with errors and warnings. An unreadable or missing
            file is reported inside that result rather than raised, so
            `beast2py validate` always prints diagnostics.
        """
        try:
            with open(xml_file, "r", encoding="utf-8") as f:
                xml_string = f.read()
        except FileNotFoundError:
            result = ValidationResult()
            result.add_error(f"XML file not found: {xml_file}")
            return result
        except OSError as exc:
            result = ValidationResult()
            result.add_error(f"Cannot read {xml_file}: {exc}")
            return result
        return XMLValidator.validate(xml_string)

    @staticmethod
    def _beast2_reported_valid(stdout: str) -> bool:
        """Decide whether BEAST2's own output affirms this file.

        A plain ``"VALID:" in output`` is satisfied by ``"INVALID: ..."``
        because VALID: is a substring of it, so the guard meant to stop an
        unrelated script exiting 0 from producing a false positive never
        actually engaged. Match line-anchored markers on stdout
        only, and treat any INVALID line as a veto.

        Args:
            stdout: The validator script's standard output.

        Returns:
            True only if a VALID line is present and no INVALID line is.
        """
        valid = False
        for line in (stdout or "").splitlines():
            stripped = line.strip()
            if stripped.startswith("INVALID:"):
                return False
            if stripped.startswith("VALID:"):
                valid = True
        return valid

    @staticmethod
    def beast2_validation_status(
        xml_file: str,
        beast2_path: str = "beast",
    ) -> Tuple[bool, str, str]:
        """Validate with BEAST2 and distinguish "failed" from "not installed".

        ``validate_with_beast2`` reports a boolean, which cannot tell a rejected
        model apart from a missing BEAST2 installation. Callers that gate on
        exit codes need that distinction: an unavailable checker must never be
        reported as a passed one.

        Args:
            xml_file: Path to the XML file.
            beast2_path: BEAST2 executable, or a beast2_validate.sh script.

        Returns:
            Tuple of (success, output, status) where status is one of
            ``"passed"``, ``"failed"`` or ``"unavailable"``.
        """
        success, output = XMLValidator.validate_with_beast2(xml_file, beast2_path=beast2_path)
        if success:
            return True, output, "passed"
        lowered = (output or "").lower()
        if (
            "setup-error:" in lowered
            or "not found" in lowered
            or "no such file" in lowered
            or "timed out" in lowered
        ):
            return False, output, "unavailable"
        return False, output, "failed"

    @staticmethod
    def _bundled_validator_script() -> Optional[Path]:
        """Locate ``beast2_validate.sh`` next to the package or in site data.

        ``pyproject.toml`` used to package only ``beast2py*``, so an installed
        distribution had no script and validation silently degraded to
        ``beast -validate`` and then to "not found".

        Returns:
            Path to the script, or None when it cannot be found.
        """
        here = Path(__file__).resolve().parent
        candidates = [
            here / "beast2_validate.sh",
            here.parent / "beast2_validate.sh",
            here.parent.parent / "beast2_validate.sh",
        ]
        # An installed distribution puts the script on PATH next to the
        # interpreter (pyproject `script-files`).
        import sys

        candidates.append(Path(sys.executable).parent / "beast2_validate.sh")
        for candidate in candidates:
            if candidate.is_file():
                return candidate
        import shutil

        on_path = shutil.which("beast2_validate.sh")
        return Path(on_path) if on_path else None

    @staticmethod
    def validate_with_beast2(
        xml_file: str,
        beast2_path: str = "beast",
    ) -> Tuple[bool, str]:
        """Validate XML using BEAST2's XMLParser (headless, no JavaFX).

        ``beast2_path`` is interpreted as follows:

        - A ``.sh``/``.bash`` file is treated as a ``beast2_validate.sh``
          script and run via ``bash``; success requires exit code 0 *and*
          the ``VALID:`` marker printed by the bundled beast2_validate.sh,
          so an unrelated script that merely exits 0 can never produce a
          false positive.
        - With the default ``"beast"`` (or a non-existent path), the bundled
          script is looked up in the project root next to the package.
        - Any other value is treated as the BEAST2 executable for a
          ``beast -validate`` fallback. For security only a fixed command
          name is ever executed: a path form must point to a file named
          ``beast`` (its directory is prepended to PATH), and a bare form
          must resolve on PATH.

        Args:
            xml_file: Path to the XML file.
            beast2_path: Path to the BEAST2 executable or a validation script.

        Returns:
            Tuple of (success, output_message).
        """
        import os
        import shutil
        import subprocess

        script: Optional[Path] = None
        beast_exec = beast2_path or "beast"

        if beast2_path and beast2_path != "beast" and os.path.isfile(beast2_path):
            if beast2_path.endswith((".sh", ".bash")):
                script = Path(beast2_path)
            # else: a BEAST2 executable path, handled by the -validate
            # fallback below.
        else:
            script = XMLValidator._bundled_validator_script()

        if script is not None:
            try:
                result = subprocess.run(
                    ["bash", str(script), xml_file],
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                output = result.stdout + result.stderr
                if result.returncode == 2:
                    # Exit code 2 is this project's documented "setup error"
                    # (no BEAST.base.jar, or no JDK 17). Mark it explicitly:
                    # a *missing* checker must never be reported as a passed
                    # -- or failed -- one, and "none was found" does not
                    # contain the substring "not found" that the status
                    # heuristic used to rely on.
                    output = f"SETUP-ERROR: {output}"
                success = (
                    result.returncode == 0
                    and XMLValidator._beast2_reported_valid(result.stdout)
                )
                return success, output
            except subprocess.TimeoutExpired:
                return False, "BEAST2 validation timed out"
            except Exception as e:
                return False, f"BEAST2 validation error: {e}"

        # Fall back to 'beast -validate'. Only the fixed command name "beast"
        # is ever executed; a custom path must be an existing executable file
        # named "beast", whose directory is prepended to PATH.
        env = dict(os.environ)
        if os.path.sep in beast_exec:
            if not (os.path.isfile(beast_exec) and os.access(beast_exec, os.X_OK)):
                return False, f"BEAST2 executable not found: {beast_exec}"
            if Path(beast_exec).name != "beast":
                return False, (
                    "Unsupported BEAST2 executable: the -validate fallback "
                    "requires a binary named 'beast' (or pass a "
                    "beast2_validate.sh script)"
                )
            env["PATH"] = str(Path(beast_exec).parent) + os.pathsep + env.get("PATH", "")
        elif beast_exec != "beast" and shutil.which(beast_exec) is None:
            return False, f"BEAST2 executable not found: {beast_exec}"
        elif shutil.which("beast") is None:
            return False, f"BEAST2 executable not found: {beast_exec}"

        try:
            result = subprocess.run(
                ["beast", "-validate", xml_file],
                capture_output=True,
                text=True,
                timeout=60,
                env=env,
            )
            success = result.returncode == 0
            output = result.stdout + result.stderr
            return success, output
        except subprocess.TimeoutExpired:
            return False, "BEAST2 validation timed out"
        except Exception as e:
            return False, f"BEAST2 validation error: {e}"

    @staticmethod
    def _collect_ids(
        elem: ET.Element,
        ids: set,
        idrefs: set,
    ) -> None:
        """Recursively collect all id and idref values."""
        if "id" in elem.attrib:
            ids.add(elem.attrib["id"])
        if "idref" in elem.attrib:
            idrefs.add(elem.attrib["idref"])

        # Also check for @ references in attribute values
        for attr_name, attr_val in elem.attrib.items():
            if attr_val.startswith("@"):
                ref_id = attr_val[1:]
                idrefs.add(ref_id)

        for child in elem:
            XMLValidator._collect_ids(child, ids, idrefs)

    @staticmethod
    def _count_ids(elem: ET.Element, id_counts: dict) -> None:
        """Count occurrences of each id."""
        if "id" in elem.attrib:
            id_val = elem.attrib["id"]
            id_counts[id_val] = id_counts.get(id_val, 0) + 1
        for child in elem:
            XMLValidator._count_ids(child, id_counts)
