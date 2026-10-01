# Security Policy

## Supported versions

| Version | Supported |
|---------|-----------|
| 0.1.x   | Yes      |
| < 0.1   | No       |

Beast2Py is pre-1.0 and has a single maintained line, 0.1.x. The current release
is 0.1.1 (tag `v0.1.1`) and nothing earlier was ever published, so there is no
older line to patch. Fixes land on `main` and in the next 0.1.x tag.

## Reporting a vulnerability

Do not open a public issue for a suspected vulnerability. It stays public for the
whole time it takes to fix the problem, and it carries crash output that may
contain paths or data from your machine.

Use GitHub's private reporting channel instead: the repository's **Security**
tab, then **Report a vulnerability**, or
<https://github.com/ZengZichao/Beast2Py/security/advisories/new>. That opens a
private Security Advisory draft only you and the maintainer can see, so the
details can be fixed and released before they are public. If you cannot reach
that form, email <zengzichao@sjtu.edu.cn> with `SECURITY` in the subject line.

## What to include

- The Beast2Py version (`beast2py --version`), Python version, operating system.
- The YAML configuration or exact command that triggers the problem. The smallest
  config that still reproduces it is worth more than the real one.
- The complete traceback, not the last line.
- Whether the BEAST2 validation gate was used (`--beast2-validate` for
  `generate` / `quick`, or `--beast2` for `validate`), and if so the BEAST2 and
  JDK versions.
- What you expected, what happened instead, and whether the issue is already
  public anywhere.

## Response time

Targets for a single-maintainer academic project, not guarantees: acknowledgement
within 7 days; a first assessment (confirmed, not reproducible, or out of scope)
within 14 days; a fix, or a written decision not to ship one, within 60 days. A
report that turns out not to be a vulnerability still gets a reply.

## What is not a vulnerability

- A divergence-time estimate you disagree with scientifically, or a calibration
  prior you think is misspecified. That is a modelling discussion: open an issue.
- A crash, traceback or rejected input on a malformed configuration. Beast2Py
  refuses unknown keys and out-of-domain priors on purpose; an uncaught exception
  is a bug, not a vulnerability. Report it on the issue tracker.
- The BEAST2 gate reporting `unavailable`, or generation stopping with exit code 2
  when BEAST2 or JDK 17 is missing. That is the documented behaviour of
  `--beast2-validate`, designed to fail loudly rather than pass.
- Anything reproducible only after editing a local checkout or patching a vendored
  class under `beast2py/tools/classes/`, and a vulnerable transitive dependency
  flagged by Dependabot. The second is still worth reporting here.

Reporting a non-vulnerability is not a mistake, and credit in the advisory is
offered if you want it. Please do not publish details of an unfixed vulnerability
before a fix or an explicit decision has been made public.
