"""Structural validation for GitHub Actions workflows.

`yaml.safe_load` is not enough: a workflow can be valid YAML and still be rejected by
GitHub, which then creates a run with zero jobs and reports the workflow by file path
instead of by name. That is exactly what happened to release.yml, where two jobs shared
the display name "publish to PyPI".

This checks the schema-level mistakes that are easy to make and invisible to a plain
YAML parse. It is a local sanity net, not a replacement for the real thing: the only
authoritative check is whether GitHub accepts the file and reports a named run.
"""

import collections
import pathlib
import sys

import yaml

TOP_LEVEL = {
    "name",
    "run-name",
    "on",
    "permissions",
    "env",
    "defaults",
    "concurrency",
    "jobs",
    "id-token",
}
JOB_KEYS = {
    "name",
    "needs",
    "if",
    "permissions",
    "environment",
    "concurrency",
    "outputs",
    "env",
    "defaults",
    "steps",
    "timeout-minutes",
    "strategy",
    "continue-on-error",
    "container",
    "services",
    "uses",
    "with",
    "secrets",
    "runs-on",
}
STEP_KEYS = {
    "id",
    "if",
    "name",
    "uses",
    "run",
    "with",
    "env",
    "continue-on-error",
    "timeout-minutes",
    "working-directory",
    "shell",
}


def validate(path: pathlib.Path) -> list[str]:
    problems: list[str] = []
    raw = path.read_text()

    if "\t" in raw:
        problems.append("contains a tab character; YAML forbids tabs for indentation")

    try:
        doc = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        return [f"not parseable as YAML: {exc}"]

    if not isinstance(doc, dict):
        return ["top level is not a mapping"]

    unknown = set(doc) - TOP_LEVEL
    # PyYAML follows YAML 1.1, where the bare key `on` is the boolean True.
    # GitHub's own parser is YAML 1.2 and reads it as the string "on", so
    # normalise rather than reporting a false positive.
    if True in doc:
        doc = {(("on" if k is True else k)): v for k, v in doc.items()}
        unknown = set(doc) - TOP_LEVEL
    if unknown:
        problems.append(f"unknown top-level keys: {sorted(unknown)}")

    if not doc.get("jobs") or not isinstance(doc["jobs"], dict):
        problems.append("no jobs defined")
        return problems

    if "on" not in doc:
        problems.append("no 'on' trigger defined")

    # The mistake that produced a zero-job run.
    display = [job.get("name", job_id) for job_id, job in doc["jobs"].items()]
    for name, count in collections.Counter(display).items():
        if count > 1:
            offenders = [j for j, v in zip(doc["jobs"], display) if v == name]
            problems.append(f"duplicate job display name {name!r} used by {offenders}")

    for job_id, job in doc["jobs"].items():
        if not isinstance(job, dict):
            problems.append(f"job {job_id!r} is not a mapping")
            continue

        unknown = set(job) - JOB_KEYS
        if unknown:
            problems.append(f"job {job_id!r} has unknown keys: {sorted(unknown)}")

        # GitHub rejects the whole workflow file if `secrets` appears in any
        # `if` expression, at job level or step level, quoted or not. Only the
        # github / needs / vars / inputs contexts are available there. The
        # symptom is distinctive and very hard to read from the run list: the
        # workflow gets registered under its file path instead of its `name`,
        # and every push produces a run containing zero jobs.
        #
        # `secrets` is fine in `env`, `run` and `with`.
        job_if = job.get("if")
        if isinstance(job_if, str) and "secrets" in job_if:
            problems.append(
                f"job {job_id!r} uses the secrets context in a job-level 'if'; "
                "use a repository variable, or move the condition to a step"
            )

        # A job must either declare its own runner or reuse a reusable workflow.
        if "uses" not in job and "runs-on" not in job:
            problems.append(f"job {job_id!r} has neither 'runs-on' nor 'uses'")

        needs = job.get("needs") or []
        needs = [needs] if isinstance(needs, str) else needs
        for dep in needs:
            if dep not in doc["jobs"]:
                problems.append(f"job {job_id!r} needs unknown job {dep!r}")

        steps = job.get("steps")
        if steps is None and "uses" not in job:
            problems.append(f"job {job_id!r} has no steps")
        for i, step in enumerate(steps or []):
            if not isinstance(step, dict):
                problems.append(f"job {job_id!r} step {i} is not a mapping")
                continue
            unknown = set(step) - STEP_KEYS
            if unknown:
                problems.append(f"job {job_id!r} step {i} has unknown keys: {sorted(unknown)}")
            if ("uses" in step) == ("run" in step):
                problems.append(f"job {job_id!r} step {i} must set exactly one of 'uses' or 'run'")
            step_if = step.get("if")
            if isinstance(step_if, str) and "secrets" in step_if:
                problems.append(
                    f"job {job_id!r} step {i} uses the secrets context in 'if'; "
                    "GitHub rejects the entire file for this, so gate on a "
                    "repository variable instead"
                )

    return problems


def main() -> int:
    targets = sorted(pathlib.Path(".github/workflows").glob("*.yml"))
    if not targets:
        print("no workflow files found", file=sys.stderr)
        return 1

    failed = False
    for path in targets:
        problems = validate(path)
        if problems:
            failed = True
            print(f"FAIL {path}")
            for p in problems:
                print(f"       - {p}")
        else:
            print(f"ok   {path}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
