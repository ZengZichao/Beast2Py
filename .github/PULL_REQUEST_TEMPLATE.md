<!-- Please fill this in: it is what the reviewer reads first. One logical change
     per pull request. A security problem does not belong in a pull request —
     see SECURITY.md. -->

## What does this change?

<!-- One or two sentences. What is different after this PR? -->

Closes #

## Why

<!-- The problem being solved, or the issue this implements; link it. If this is
     science-facing (a new model, a different XML shape, a changed validation
     rule), say what the previous behaviour did and why it was wrong. -->

## Checklist

- [ ] `flake8 beast2py tests --max-line-length=100` is clean, and
      `pytest tests/ -v -rs` passes locally.
- [ ] A test covers the change, and a test that would have caught the old
      behaviour now fails without this change.
- [ ] If `beast2py/tools/Beast2Validator.java` changed, the recompiled
      `beast2py/tools/classes/Beast2Validator.class` is committed in this PR too.
      The committed bytecode is what CI and the wheel actually run; stale
      bytecode makes the change silently do nothing.
- [ ] If a documented message, command, flag, path or count changed, that is
      reflected in `README.md`, `docs/` and `tests/test_release_integrity.py`.
- [ ] Every `docs/*_en.md` change is mirrored in the matching `docs/*_zh.md`, and
      `README.md` and `README_zh.md` are changed together.
- [ ] `CHANGELOG.md` has an entry under `## [Unreleased]` for anything a user
      would notice.
- [ ] No new third-party binary, jar or `.class` file, and no GPL-only or
      LGPL-incompatible code. If you need one, say so here: that needs an issue
      first, not a PR.

## Verification

Paste the summary line from `pytest tests/ -rs` (or `pytest tests/ -v -rs`):

```
N passed, M skipped in Xs
```

- Did the BEAST2 integration tests **execute or skip** on your machine? State
  which, and paste the `-rs` skip reasons if any:

  - [ ] executed (BEAST2 2.7.x + JDK 17 present)
  - [ ] skipped (no BEAST2 / no JDK 17 — nothing is verified about the generated
        XML in that case)

- If the change touches XML generation, did you run the third gate yourself with
  `beast2py generate --config <your config> --output out.xml --beast2-validate`?
- If a generated XML, fingerprint or methods file changed, name the example you
  regenerated and say whether `examples/output_basic.xml` and its
  `.fingerprint.json` / `.methods.tex` sidecars are still in sync.
- A green CI run does not prove the BDSKY or nested-sampling cases ran: a failed
  add-on download downgrades them to skips and the job stays green.
