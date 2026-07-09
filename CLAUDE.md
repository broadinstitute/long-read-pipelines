# Working in this repository

## Start here

Before doing anything else, read [`docs/development_guide/README.md`](docs/development_guide/README.md)
and every file it links to:

- [`docs/development_guide/repo_structure.md`](docs/development_guide/repo_structure.md) — where things live
- [`docs/development_guide/WDL_STYLE_RULES.md`](docs/development_guide/WDL_STYLE_RULES.md) — authoritative WDL style spec
- [`docs/development_guide/docker_style_guide.md`](docs/development_guide/docker_style_guide.md) — Docker conventions

Those documents are the source of truth for conventions. The notes below are practical
tips accumulated while working in the repo; they supplement, not replace, the guide.

## Validating WDL

Validate every changed `.wdl` with **both** tools before committing (see the
"WDL Validation" section of the development guide):

- `miniwdl check --strict <file.wdl>` — fast, and lints embedded `command` blocks via shellcheck.
- `java -jar womtool.jar validate <file.wdl>` — the Terra/Cromwell source of truth. Download
  `womtool-<version>.jar` from the Cromwell GitHub releases.

miniwdl is more lenient than womtool: it accepts some things Cromwell rejects (e.g. a
parenthesized boolean inside a `~{true=... false=... expr}` placeholder — hoist such an
expression into a `Boolean` declaration and interpolate the bare variable). A WDL is not
"valid" until womtool passes. Validate the top-level workflow so imported task libraries are
checked transitively.

## WDL conventions (quick reference; the style rules are authoritative)

- Structs import with **no alias**: `import "../../structs/Structs.wdl"` (path depth varies).
- `command <<< ... >>>` heredoc, first line `set -euxo pipefail`.
- Runtime via the `RuntimeAttr` struct + `runtime_attr_override` pattern, with defaults
  `boot_disk_gb: 25`, `preemptible_tries: 1` (2 for cheap/idempotent tasks), `max_retries: 1`.
- Prefer keeping `task` blocks in a task-library `.wdl` (e.g. `wdl/tasks/Utility/…`) separate
  from the `workflow` file (`wdl/pipelines/…`).
- Register every **new pipeline workflow** in `.dockstore.yml` (task libraries are not registered).

## Streaming from GCS / images

- To stream a `gs://` file into bcftools/htslib without localizing: mark the `File` input
  `localization_optional: true` in `parameter_meta`, then export a token in the `command`.
  On an image with `gcloud` (e.g. `lr-basic`): `GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`
  (use the `application-default` variant — plain `gcloud auth print-access-token` has no active
  account on a Cromwell VM). On images without gcloud, hit the metadata server directly:
  `GCS_OAUTH_TOKEN=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token | grep -o '"access_token":"[^"]*"' | cut -d'"' -f4)`.
- Image notes: `us.gcr.io/broad-dsp-lrma/lr-basic:0.1.x` has bcftools/samtools (libcurl-enabled),
  `curl`, and `gcloud`, but **no python3** — use `python:3.11-slim` for pure-stdlib python tasks.
  `firecloud.api` (FISS) is only in `us.gcr.io/broad-dsp-lrma/lr-backup-workspace:0.0.1`.
- Embedded `python3 <<CODE ... CODE` heredocs work because Cromwell dedents the whole command
  block; keep the python body's indentation internally consistent.

## Terra data tables

- Auto-detect namespace/workspace inside a task via `WORKSPACE_NAMESPACE`/`WORKSPACE_NAME`
  env vars, else match `GOOGLE_PROJECT` (metadata server) against `fapi.list_workspaces`.
- Entity references (clickable links) are expressible in a flexible-import TSV as a JSON cell:
  `{"entityType":"<table>","entityName":"<id>"}` (arrays via `{"itemsType":"EntityReference","items":[...]}`).

## Git workflow

- Branch names start with the author's initials, e.g. `jts_<short_description>`.
- Commit after every material change. Never use `--no-verify`.
- The Git-LFS `post-commit` warning ("git-lfs not found") is benign; the commit still succeeds.
