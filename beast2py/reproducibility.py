"""Reproducibility framework for Beast2Py.

Core innovation module containing:
- FingerprintGenerator: Generates unique analysis identifiers
- MethodsGenerator: Automatically generates LaTeX methods descriptions
- PipelineGenerator: Generates Snakemake/Nextflow pipeline files
"""

from __future__ import annotations

import math
import hashlib
import re
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from . import __version__
from .models import (
    BEASTConfig,
    ClockModelType,
    TreePriorType,
)


class FingerprintGenerator:
    """Generate analysis fingerprints for reproducibility.

    The identifier is a pure function of the *analysis*: configuration content
    plus the content of the alignments. It deliberately contains no calendar
    date, no time zone and not the output file name, so the same YAML plus the
    same data yields the same identifier whenever and wherever it is run, while
    a changed model component, calibration, MCMC setting *or sequence* changes
    it ( and the's.
    """

    # 12 hex characters = 48 bits of digest. By the birthday bound, a 50%
    # chance of one collision needs ~2^24 (≈1.7e7) distinct analyses, so the
    # collision risk is irrelevant at the scale of published phylogenies.
    HASH_HEX_CHARS = 12

    @classmethod
    def collision_bits(cls) -> int:
        """Return the number of digest bits the identifier carries."""
        return 4 * cls.HASH_HEX_CHARS

    @classmethod
    def collision_bound_50pct(cls) -> float:
        r"""Return n at which a truncated digest has a 50% chance of one collision.

        For n values hashed uniformly into ``2**bits`` slots,
        ``P(>=1 collision) = 1 - exp(-n(n - 1) / 2**(bits + 1))``, which reaches
        0.5 at ``n = sqrt(2 * 2**bits * ln 2)`` -- that is 1.1774 * 2**(bits/2),
        about 2.0e7 for the 48 bits emitted here.  ``sqrt(2**48)`` on its own is
        only a 39% chance, so the two must not be conflated; the earlier wording
        in this file did exactly that.
        """
        return math.sqrt(2 * 2 ** cls.collision_bits() * math.log(2))

    @classmethod
    def collision_note(cls) -> str:
        """Return the collision caveat written into the sidecar."""
        bits = cls.collision_bits()
        return (
            f"The identifier truncates the configuration digest to {bits} bits; "
            f"a 50% chance of one collision needs sqrt(2 * 2**{bits} * ln 2), "
            f"about {cls.collision_bound_50pct():.2e} distinct analyses, so it is "
            "a change-detection token rather than a globally unique key."
        )

    @staticmethod
    def _canonical(obj):
        """Serialise anything ``json`` cannot, without inheriting container order.

        ``default=str`` was the fallback here and is not order-stable:
        ``str({1, 2})`` differs from ``str({2, 1})``, so a set ever reaching
        ``to_dict()`` would make the digest depend on ``PYTHONHASHSEED`` and
        break the byte-identical-rerun claim. Sorting removes that hazard.
        """
        if isinstance(obj, (set, frozenset)):
            return sorted(str(x) for x in obj)
        if isinstance(obj, dict):
            return {str(k): obj[k] for k in sorted(obj, key=str)}
        return str(obj)

    @staticmethod
    def config_hash(config: BEASTConfig) -> str:
        """Return the first 12 hex digits of the configuration digest.

        Args:
            config: BEAST configuration.

        Returns:
            12-character lowercase hex string.
        """
        config_str = json.dumps(
            config.to_dict(),
            sort_keys=True,
            default=FingerprintGenerator._canonical,
            ensure_ascii=False,
        )
        return hashlib.sha256(config_str.encode("utf-8")).hexdigest()[
            : FingerprintGenerator.HASH_HEX_CHARS
        ]

    @staticmethod
    def data_hash(config: BEASTConfig) -> str:
        """Return a 16-hex digest over the content of every alignment.

        The digest covers alignment *content*, not just the file paths, so
        replacing a FASTA's sequences while keeping its name changes the
        analysis identity instead of leaving it untouched.

        Args:
            config: BEAST configuration.

        Returns:
            First 16 hex digits of the digest over the per-alignment digests.
        """
        hasher = hashlib.sha256()
        seen = set()
        for part in sorted(config.partitions, key=lambda p: p.id):
            if part.alignment.id in seen:
                continue
            seen.add(part.alignment.id)
            hasher.update(part.alignment.id.encode("utf-8"))
            hasher.update(b"=")
            hasher.update(part.alignment.content_digest().encode("utf-8"))
            hasher.update(b";")
        return hasher.hexdigest()[:16]

    @staticmethod
    def generate_fingerprint(config: BEASTConfig) -> str:
        """Generate the unique analysis identifier.

        Format: ``B2P-{hash12}-{version}``

        Args:
            config: BEAST configuration.

        Returns:
            Fingerprint string, e.g. ``B2P-a3f7b2c9d1e8-0.1.0``.
        """
        tool_version = config.metadata.get("tool_version") or __version__
        return f"B2P-{FingerprintGenerator.config_hash(config)}-{tool_version}"

    @staticmethod
    def generate_comment(config: BEASTConfig) -> str:
        """Build the XML comment line embedded in generated files.

        Contains the deterministic identifier and the data digest, but no
        generation timestamp: a timestamp inside the artefact would make the
        "reproducible product" differ byte-wise on every rerun, contradicting
        the claim it is meant to document. The wall-clock time is recorded in
        the JSON sidecar instead.

        Args:
            config: BEAST configuration.

        Returns:
            Comment body, without the surrounding ``<!-- -->``.
        """
        return (
            f" Analysis Fingerprint: {FingerprintGenerator.generate_fingerprint(config)}"
            f" | Data: {FingerprintGenerator.data_hash(config)} "
        )

    @staticmethod
    def generate_fingerprint_dict(
        config: BEASTConfig, xml_content: Optional[str] = None
    ) -> Dict[str, Any]:
        """Generate a detailed fingerprint dictionary.

        Args:
            config: BEAST configuration.
            xml_content: The generated XML, when available. Its digest is added
                so the sidecar certifies the artefact BEAST2 actually runs, not
                just the YAML it came from.

        Returns:
            Dictionary with fingerprint details, including the run time (kept
            out of the identifier) and the alignment content digest.
        """
        fingerprint = FingerprintGenerator.generate_fingerprint(config)
        config_str = json.dumps(
            config.to_dict(),
            sort_keys=True,
            default=FingerprintGenerator._canonical,
            ensure_ascii=False,
        )
        full_hash = hashlib.sha256(config_str.encode("utf-8")).hexdigest()

        result = {
            "fingerprint": fingerprint,
            "config_hash": full_hash,
            "full_hash": full_hash,
            "data_hash": FingerprintGenerator.data_hash(config),
            "hash_bits": 4 * FingerprintGenerator.HASH_HEX_CHARS,
            "collision_note": FingerprintGenerator.collision_note(),
            "tool_version": config.metadata.get("tool_version") or __version__,
            "beast2_version": config.metadata.get("beast2_version", "2.7.8"),
            "generation_time": datetime.now().isoformat(),
            "analysis_name": config.metadata.get("analysis_name", ""),
        }
        if xml_content is not None:
            # The sidecar used to certify only the YAML, so a hand-edited XML --
            # or a regression in the writer -- still "matched" its fingerprint.
            # What BEAST2 actually runs is the XML.
            result["xml_digest"] = hashlib.sha256(
                xml_content.encode("utf-8")
            ).hexdigest()
        return result

    @staticmethod
    def xml_digest(xml_content: str) -> str:
        """SHA-256 of generated XML text.

        Args:
            xml_content: The XML as written to disk.

        Returns:
            Hex digest of the exact bytes the file contains.
        """
        return hashlib.sha256(xml_content.encode("utf-8")).hexdigest()

    @staticmethod
    def write_fingerprint_file(
        config: BEASTConfig, output_file: str, xml_content: Optional[str] = None
    ) -> str:
        """Write fingerprint information to a JSON file.

        Args:
            config: BEAST configuration.
            output_file: Output file path.
            xml_content: The generated XML, when available.

        Returns:
            Path to the written file.
        """
        fp_dict = FingerprintGenerator.generate_fingerprint_dict(config, xml_content)
        Path(output_file).write_text(json.dumps(fp_dict, indent=2), encoding="utf-8")
        return output_file


class MethodsGenerator:
    """Automatically generate LaTeX methods descriptions."""

    # Reference citations
    CITATIONS = {
        "beast2": "Bouckaert et al. (2014) PLoS Comput Biol 10: e1003537",
        "hky": "Hasegawa et al. (1985) J Mol Evol 22: 160-174",
        "gtr": "Tavare (1986) Lect Math Life Sci 17: 57-86",
        "jc69": "Jukes & Cantor (1969) In: Mammalian Protein Metabolism, pp. 21-132",
        "tn93": "Tamura & Nei (1993) Mol Biol Evol 10: 512-526",
        "yule": "Yule (1925) Philos Trans R Soc Lond B 213: 21-87",
        "birth_death": "Gernhard (2008) J Theor Biol 253: 769-778",
        "coalescent": "Kingman (1982) Stoch Process Their Appl 13: 235-248",
        "ucln": "Drummond et al. (2006) PLoS Biol 4: e88",
        "rlc": "Drummond & Suchard (2010) BMC Biol 8: 114",
        "mrca": "Heled & Drummond (2012) Syst Biol 61: 138-149",
        "calibrated_yule": "Heled & Drummond (2012) Syst Biol 61: 138-149",
        "bsp": "Drummond et al. (2005) Mol Biol Evol 22: 1185-1192",
        "ebsp": "Heled & Drummond (2008) BMC Evol Biol 8: 289",
        "gamma": "Yang (1994) J Mol Evol 39: 306-314",
        "calibration_best_practice": "Rieux & Balloux (2016) Mol Ecol 25: 1911-1924",
        "reproducible": "Sandve et al. (2013) PLoS Comput Biol 9: e1003285",
    }

    @staticmethod
    def cite(key: str) -> str:
        """In-text form of a reference, e.g. "Bouckaert et al. 2014".

        The paragraph is pasted into a manuscript, so citations must follow the
        journal's author-year form; the full strings are emitted after it.
        """
        full = MethodsGenerator.CITATIONS.get(key, "")
        m = re.match(r"^(.*?) \((\d{4})\)", full)
        if not m:
            return full
        return m.group(1).replace(" & ", " and ") + " " + m.group(2)

    @staticmethod
    def references(paragraph: str) -> str:
        """The full entries behind the short citations used in ``paragraph``."""
        seen, lines = set(), []
        for key, full in MethodsGenerator.CITATIONS.items():
            short = MethodsGenerator.cite(key)
            if short in paragraph and full not in seen:
                seen.add(full)
                lines.append("% " + full.replace(" & ", " and "))
        if not lines:
            return ""
        return ("\n\n% Cited author-year keys. The reference list must name every author and"
                "\n% give the full title, so expand each line below before submitting:\n"
                + "\n".join(lines))

    @staticmethod
    def generate_methods(config: BEASTConfig) -> str:
        """Generate a LaTeX-formatted methods description.

        Args:
            config: BEAST configuration.

        Returns:
            LaTeX text string.
        """
        # Analysis overview
        beast2_ver = config.metadata.get("beast2_version", "2.7.8")

        # Build substitution model description
        subst_desc = MethodsGenerator._describe_substitution_models(config)

        # Build clock model description
        clock_desc = MethodsGenerator._describe_clock_models(config)

        # Build tree prior description
        tree_prior_desc = MethodsGenerator._describe_tree_prior(config)

        # Build calibration description
        cal_desc = MethodsGenerator._describe_calibrations(config)

        # Build MCMC description
        mcmc_desc = MethodsGenerator._describe_mcmc(config)

        fingerprint = FingerprintGenerator.generate_fingerprint(config)

        # Combine into full methods paragraph
        methods_text = (
            f"Divergence times were estimated using BEAST2 v{beast2_ver} "
            f"({MethodsGenerator.cite('beast2')}). "
            f"{subst_desc} "
            f"{clock_desc} "
            f"{tree_prior_desc} "
            f"{cal_desc} "
            f"{mcmc_desc} "
            f"The analysis was configured using Beast2Py v{__version__} "
            f"({MethodsGenerator.cite('reproducible')}; "
            f"Analysis Fingerprint: {fingerprint})."
        )

        return methods_text + MethodsGenerator.references(methods_text)

    @staticmethod
    def _describe_substitution_models(config: BEASTConfig) -> str:
        """Describe the substitution models used."""
        parts: List[str] = []

        for p in config.partitions:
            sm_type = p.site_model.substitution_model.get("type", "hky")
            sm_type_lower = sm_type.lower()

            # Model name mapping
            model_names = {
                "hky": "HKY",
                "gtr": "GTR",
                "jc69": "Jukes-Cantor (JC69)",
                "tn93": "TN93",
                "sym": "SYM",
                "tim": "TIM",
                "tvm": "TVM",
                "wag": "WAG",
                "jtt": "JTT",
                "dayhoff": "Dayhoff",
                "blosum62": "BLOSUM62",
            }
            model_name = model_names.get(sm_type_lower, sm_type.upper())
            citation = MethodsGenerator.cite(sm_type_lower)

            gamma_cats = p.site_model.gamma_categories
            gamma_desc = ""
            if gamma_cats > 0:
                gamma_desc = (
                    f" with gamma-distributed rate heterogeneity "
                    f"({gamma_cats} categories; {MethodsGenerator.cite('gamma')})"
                )

            prop_inv = p.site_model.proportion_invariant
            inv_desc = ""
            if prop_inv and float(prop_inv.value) > 0:
                inv_desc = " and a proportion of invariant sites"

            citation_str = f" ({citation})" if citation else ""
            parts.append(
                f"Sequence evolution for partition '{p.id}' was modeled using a "
                f"{model_name} substitution model{citation_str}{gamma_desc}{inv_desc}"
            )

        return ". ".join(parts) + "."

    @staticmethod
    def _describe_clock_models(config: BEASTConfig) -> str:
        """Describe the clock models used."""
        clock_descs: Dict[str, str] = {}

        for p in config.partitions:
            if p.clock_model.linked_to:
                continue

            ct = p.clock_model.type
            if ct == ClockModelType.STRICT:
                clock_descs[p.id] = "a strict molecular clock"
            elif ct == ClockModelType.UCLN:
                clock_descs[p.id] = (
                    f"an uncorrelated relaxed clock with lognormally distributed rates "
                    f"({MethodsGenerator.cite('ucln')})"
                )
            elif ct == ClockModelType.UCE:
                clock_descs[p.id] = (
                    f"an uncorrelated relaxed clock with exponentially distributed rates "
                    f"({MethodsGenerator.cite('ucln')})"
                )
            elif ct == ClockModelType.RLC:
                clock_descs[p.id] = (
                    f"a random local clock model " f"({MethodsGenerator.cite('rlc')})"
                )

        if not clock_descs:
            return ""
        parts = []
        for pid, desc in clock_descs.items():
            if len(clock_descs) == 1:
                parts.append(desc[0].upper() + desc[1:] + " was assumed")
            else:
                parts.append(f"for partition '{pid}', {desc} was used")

        if len(parts) == 1:
            return parts[0] + "."
        return ". ".join(parts) + "."

    @staticmethod
    def _describe_tree_prior(config: BEASTConfig) -> str:
        """Describe the tree prior."""
        tp = config.tree_prior_type

        if tp == TreePriorType.YULE:
            return (
                f"A Yule speciation process ({MethodsGenerator.cite('yule')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.BIRTH_DEATH:
            return (
                f"A birth-death speciation process "
                f"({MethodsGenerator.cite('birth_death')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.COALESCENT_CONSTANT:
            return (
                f"A coalescent process with constant population size "
                f"({MethodsGenerator.cite('coalescent')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.COALESCENT_EXPONENTIAL:
            return (
                f"A coalescent process with exponential population growth "
                f"({MethodsGenerator.cite('coalescent')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.BAYESIAN_SKYLINE:
            return (
                f"A Bayesian skyline coalescent model "
                f"({MethodsGenerator.cite('bsp')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.CALIBRATED_YULE:
            return (
                f"A calibrated Yule model "
                f"({MethodsGenerator.cite('calibrated_yule')}) "
                f"was used as the tree prior, integrating calibration information "
                f"directly into the tree prior."
            )
        elif tp == TreePriorType.EBSP:
            return (
                f"An extended Bayesian skyline plot "
                f"({MethodsGenerator.cite('ebsp')}) "
                f"was used as the tree prior."
            )
        elif tp == TreePriorType.BD_SKYLINE_SERIAL:
            return (
                "A birth-death skyline tree prior for serially sampled data "
                "(bdsky add-on) was used as the tree prior."
            )
        else:
            return f"A {tp.value.replace('_', ' ')} tree prior was used."

    @staticmethod
    def _describe_calibrations(config: BEASTConfig) -> str:
        """Describe the calibration points."""
        if not config.calibrations:
            return "No calibration points were applied."

        parts: List[str] = []
        n_cals = len(config.calibrations)

        parts.append(
            f"{n_cals} calibration point{'s were' if n_cals > 1 else ' was'} applied "
            f"({MethodsGenerator.cite('mrca')}; "
            f"{MethodsGenerator.cite('calibration_best_practice')})"
        )

        for cal in config.calibrations:
            taxa_str = ", ".join(cal.taxa) if cal.taxa else "all taxa (root)"

            dist_desc = "an unspecified distribution"
            if cal.distribution:
                dt = cal.distribution.type.lower()
                params = cal.distribution.parameters

                if dt == "normal":
                    dist_desc = (
                        f"a Normal distribution (mean={params.get('mean', '?')}, "
                        f"sigma={params.get('sigma', '?')})"
                    )
                elif dt in ("lognormal", "log_normal"):
                    dist_desc = (
                        f"a LogNormal distribution (M={params.get('M', '?')}, "
                        f"S={params.get('S', '?')})"
                    )
                elif dt == "uniform":
                    dist_desc = (
                        f"a Uniform distribution (lower={params.get('lower', '?')}, "
                        f"upper={params.get('upper', '?')})"
                    )
                elif dt == "exponential":
                    dist_desc = f"an Exponential distribution (mean={params.get('mean', '?')})"
                elif dt == "gamma":
                    dist_desc = (
                        f"a Gamma distribution (alpha={params.get('alpha', '?')}, "
                        f"beta={params.get('beta', '?')})"
                    )

                if cal.distribution.offset and cal.distribution.offset > 0:
                    dist_desc += f" with offset {cal.distribution.offset}"

            ref_str = ""
            if cal.provenance and cal.provenance.reference:
                ref_str = f" based on {cal.provenance.reference}"

            cal_type_str = ""
            if cal.provenance and cal.provenance.calibration_type:
                cal_type_str = f" ({cal.provenance.calibration_type} boundary)"

            parts.append(
                f"the MRCA of {taxa_str} was calibrated with {dist_desc}{ref_str}{cal_type_str}"
            )

        return ": ".join([parts[0]] + [f"{p}" for p in parts[1:]]) + "."

    @staticmethod
    def _describe_mcmc(config: BEASTConfig) -> str:
        """Describe the MCMC settings."""
        chain_length = config.mcmc.chain_length
        log_every = config.loggers.trace_log.get("log_every", 1000)

        text = (
            f"MCMC chains were run for {chain_length:,} generations, "
            f"sampling every {log_every:,} generations"
        )

        if config.mcmc.mcmc_type.value == "nested_sampling":
            text += (
                f" using nested sampling (particle count={config.mcmc.particle_count}, "
                f"sub-chain length={config.mcmc.sub_chain_length})"
            )

        return text + "."


def _relative_to(path: str, base: str) -> str:
    """``path`` as seen from ``base``, unchanged when no relative form exists."""
    try:
        return os.path.relpath(path, base)
    except ValueError:  # different drive on Windows
        return path


class PipelineGenerator:
    """Generate Snakemake/Nextflow pipeline files."""

    @staticmethod
    def generate_snakemake(
        config: BEASTConfig,
        output_dir: str,
        config_file: str = "config.yaml",
    ) -> str:
        """Generate a Snakemake pipeline file.

        Args:
            config: BEAST configuration.
            output_dir: Output directory.
            config_file: Path to the YAML config file.

        Returns:
            Path to the generated Snakefile.
        """
        os.makedirs(output_dir, exist_ok=True)
        analysis_name = config.metadata.get("analysis_name", "analysis")
        output_xml = config.metadata.get("output_file", f"{analysis_name}.xml")

        # Track the actual sequence files so Snakemake re-generates the XML
        # when they change (they may live relative to the config file).  Paths
        # are written as seen from the Snakefile, because Snakemake resolves
        # inputs against the workflow directory; an absolute path would record
        # the host that generated the pipeline and break everywhere else.
        seq_files = sorted({a.source_file for a in config.all_alignments if a.source_file})
        if seq_files:
            seq_list = ",\n            ".join(
                repr(_relative_to(str(f), str(output_dir))) for f in seq_files)
            seq_input = ",\n        sequences = [\n            " + seq_list + "\n        ]"
        else:
            seq_input = ""

        # BEAST2 chooses the names of its own outputs from the XML's loggers, and
        # the $(seed) macro inside those names is expanded by the launcher at run
        # time.  The rule therefore claims a sentinel it creates itself and states
        # the real names in a comment instead of inventing predictable ones.
        trace_name = str(config.loggers.trace_log.get("file_name", ""))
        tree_name = str(config.loggers.tree_log.get("file_name", ""))

        snakefile_content = f"""# Snakemake pipeline for Beast2Py
# Analysis: {analysis_name}
# Generated by Beast2Py v{__version__}

configfile: "{config_file}"

# Rule 1: Generate BEAST2 XML
rule generate_xml:
    output:
        "{output_xml}"
    input:
        config_file = "{config_file}"{seq_input}
    shell:
        "beast2py generate --config {{input.config_file}} --output {{output}}"

# Rule 2: Run BEAST2
# The run writes the logger files named in the XML -- trace "{trace_name}",
# trees "{tree_name}" -- where $(seed) is substituted by the beast launcher,
# so this rule publishes a completion marker rather than guessing those names.
rule run_beast:
    output:
        "beast2_run.done"
    input:
        xml = "{output_xml}"
    shell:
        "beast -overwrite -threads 4 {{input.xml}} && touch {{output}}"

# Rule 3: Run calibration diagnostics
rule diagnose:
    output:
        report = "diagnostic_report.html"
    input:
        config_file = "{config_file}"
    shell:
        "beast2py diagnose --config {{input.config_file}} "
        "--report {{output.report}} --sensitivity --visualize"

# Rule 4: Generate methods description
rule methods:
    output:
        tex = "methods.tex"
    input:
        config_file = "{config_file}"
    shell:
        "beast2py methods --config {{input.config_file}} --output {{output.tex}}"

# Rule: All
rule all:
    input:
        rules.generate_xml.output,
        rules.run_beast.output,
        rules.diagnose.output,
        rules.methods.output
"""

        snakefile_path = os.path.join(output_dir, "Snakefile")
        Path(snakefile_path).write_text(snakefile_content, encoding="utf-8")

        return snakefile_path

    @staticmethod
    def generate_nextflow(
        config: BEASTConfig,
        output_dir: str,
        config_file: str = "config.yaml",
    ) -> str:
        """Generate a Nextflow pipeline file.

        Args:
            config: BEAST configuration.
            output_dir: Output directory.
            config_file: Path to the YAML config file.

        Returns:
            Path to the generated Nextflow file.
        """
        os.makedirs(output_dir, exist_ok=True)
        analysis_name = config.metadata.get("analysis_name", "analysis")
        output_xml = config.metadata.get("output_file", f"{analysis_name}.xml")
        trace_name = str(config.loggers.trace_log.get("file_name", ""))
        tree_name = str(config.loggers.tree_log.get("file_name", ""))

        nextflow_content = f'''// Nextflow pipeline for Beast2Py
// Analysis: {analysis_name}
// Generated by Beast2Py v{__version__}

params.config_file = "{config_file}"
params.output_xml = "{output_xml}"

process generate_xml {{
    output path "{output_xml}"

    script:
    """
    beast2py generate \\
        --config ${{params.config_file}} \\
        --output ${{params.output_xml}}
    """
}}

// The run writes the logger files named in the XML -- trace "{trace_name}",
// trees "{tree_name}" -- where $(seed) is substituted by the beast launcher,
// so the process publishes a completion marker rather than guessing those names.
process run_beast {{
    input path xml
    output path "beast2_run.done"

    script:
    """
    beast -overwrite -threads 4 ${{xml}}
    touch beast2_run.done
    """
}}

process diagnose {{
    output path "diagnostic_report.html"

    script:
    """
    beast2py diagnose \\
        --config ${{params.config_file}} \\
        --report diagnostic_report.html \\
        --sensitivity --visualize
    """
}}

process methods {{
    output path "methods.tex"

    script:
    """
    beast2py methods \\
        --config ${{params.config_file}} \\
        --output methods.tex
    """
}}

workflow {{
    generate_xml()
    run_beast(generate_xml.out)
    diagnose()
    methods()
}}
'''

        nextflow_path = os.path.join(output_dir, "main.nf")
        Path(nextflow_path).write_text(nextflow_content, encoding="utf-8")

        return nextflow_path
