# WDL Style Rules

Authoritative ruleset for writing and reformatting WDL files in this repo, derived from the conventions in `broadinstitute/long-read-pipelines`. These rules are intended to be machine-readable: an LLM (or a human reviewer) following them should produce WDL files that are stylistically indistinguishable from the canonical examples.

---

## 1. Repo layout

```
wdl/
├── structs/Structs.wdl                       # shared RuntimeAttr struct, lives at one canonical path
├── tasks/<Category>/<Name>.wdl               # task modules grouped by function (Alignment, QC, Utility, ...)
└── pipelines/<Tech>/<Category>/<Name>.wdl    # workflows; <Tech> = PacBio | ONT | ILMN | TechAgnostic
```

Filename = task or workflow name, CamelCase, `.wdl`.

Keep task blocks in a separate `.wdl` file from the workflow block that calls them (tasks under `tasks/`, workflows under `pipelines/`). This keeps workflows modular, readable, and maintainable, and lets multiple workflows share one task module.

---

## 2. Boilerplate (every file)

```wdl
version 1.0

import "../../structs/Structs.wdl"                        # task file (depth 2)
# or
import "../../../tasks/Utility/Utils.wdl" as Utils        # pipeline file (depth 3), always aliased
```

- `version 1.0` always.
- Tasks import `Structs.wdl` only.
- Pipelines import every task module they call, **with an `as <Alias>`**. Common aliases: `PB`, `Utils`, `FF` (Finalize), `AM` (AlignedMetrics), etc.
- Relative paths only. Count `../` from the file's depth.

---

## 3. Task skeleton (canonical order)

```wdl
task ExampleTask {

    meta {
        description: "One-line summary of what task does."

        # Required when task wraps a specific named tool:
        tool:          "samtools sort"
        tool_version:  "1.23.1"
        tool_url:      "https://www.htslib.org/"
        tool_citation: "Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. GigaScience. 2021;10(2):giab008."

        # Optional:
        author: "Jonn Smith"
        email:  "jonn@broadinstitute.org"

        # Required when task has an output {} block — see §16:
        outputs: {
            output_bam: "Coordinate-sorted BAM",
            output_bai: "BAM index (.bam.bai)"
        }
    }

    parameter_meta {
        input_bam:             "BAM to operate on"
        prefix:                "basename for outputs"
        extra_args:            "Additional command-line args appended verbatim to the samtools sort invocation"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam
        String prefix

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(5.0 * size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        # ---- Resource detection (required preamble) ----
        NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
        RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

        # Reserve 1 GB for OS + container overhead.
        USABLE_RAM_GB=$((RAM_IN_GB - 1))
        [[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

        # Per-thread RAM (used by samtools sort -m, etc.)
        MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
        [[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

        # Java heap (used by GATK/Picard). Same as USABLE_RAM_GB; alias for clarity.
        JAVA_MEM_GB=${USABLE_RAM_GB}

        echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
        # ---- end preamble ----

        samtools sort -@ ${NUM_CPUS} -m ${MEM_PER_THREAD_GB}G ~{extra_args} ~{input_bam} -o ~{prefix}.bam
        samtools index -@ ${NUM_CPUS} ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
```

### Section order inside `task`

1. `meta { ... }` (description; optional tool citation; optional author/email; required `outputs:` when task has outputs)
2. `parameter_meta { ... }` — entry per input, including `extra_args` and `runtime_attr_override`
3. `input { ... }` — required first, optional (`Type?`) after a blank line, `String extra_args = ""` then `RuntimeAttr? runtime_attr_override` last
4. Private declarations (e.g. `Int disk_size = ...`)
5. `command <<< ... >>>` — first line `set -euxo pipefail`, then the resource-detection preamble (§14)
6. `output { ... }`
7. `#########################` separator (literal, exactly as shown)
8. `RuntimeAttr default_attr = object { ... }`
9. `RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])`
10. `runtime { ... }` — column-aligned `select_first` per field

### Tool citation rules (in `meta`)

- **Mandatory** when task wraps a specific named tool (samtools, minimap2, GATK, bcftools, fastqc, etc.). Skip only for pure-utility tasks (`ChunkManifest`, `MakeChrIntervalList`, ad-hoc bash glue).
- Keys:
  - `tool` — canonical name of primary binary invoked (matches `~{extra_args}` target from §15).
  - `tool_version` — pinned version. Must match the docker tag.
  - `tool_url` — upstream homepage or GitHub.
  - `tool_citation` — full bibliographic citation if the tool has a publication. Drop if none.
- When task chains multiple tools (sort + index), cite the **primary** tool only — same one `extra_args` flows to.
- Update `tool_version` whenever the docker tag bumps.

---

## 4. `parameter_meta` rules

- Document every input. Either short string `name: "desc"` or object form for tagged inputs:
  ```wdl
  file: {
      description: "file to finalize",
      localization_optional: true
  }
  ```
- Align colons vertically when reasonable (cosmetic — not enforced).
- `extra_args` documented as `"Additional command-line args appended verbatim to the <tool> invocation"` — replace `<tool>` with the actual binary name.
- `runtime_attr_override` always documented as `"Override the default runtime attributes"`.

---

## 5. `command` block rules

- `command <<< ... >>>` (heredoc form, not `command { ... }`). Required for `~{var}` substitution.
- First line: `set -euxo pipefail`. Immediately followed by the §14 resource-detection preamble.
- Interpolate WDL values with `~{var}`, never `${var}` (`${}` reserved for bash inside the heredoc).
- For Python multi-line: use `python3 <<CODE ... CODE` heredoc; `~{var}` still expands inside.

---

## 6. Disk sizing convention

**Default: 5× total input file size, plus a 10 GB constant for output + scratch.**

```wdl
Int disk_size = 10 + ceil(5.0 * size(input_bam, "GB"))
```

Multiple inputs:

```wdl
Int disk_size = 10 + ceil(5.0 * (size(input_bam, "GB") + size(ref_fasta, "GB")))
```

Array inputs:

```wdl
Int disk_size = 10 + ceil(5.0 * size(reads, "GB"))
```

### Rules

- **Constant** `10` GB floor covers output + tmp + docker layer.
- **Multiplier 5.0** = input localization + intermediates + output + scratch slack. Default for any new task.
- **Override the 5× when task obviously breaks it:**
  - Alignment / sort / BAM-rewriting with large intermediate SAMs: bump to ~10× (`1 + 10*ceil(...)`).
  - Assembly / k-mer counting: hard-code (e.g. 3000 GB SSD).
  - Pure-passthrough (rename, index-only): fixed small disk (`disk_gb: 10`).
- Always `ceil(...)`.
- `5.0 * size(...)` returns `Float`; `ceil` converts to `Int`.
- Plug into `default_attr`:
  ```wdl
  disk_gb: disk_size,
  ```

---

## 7. Runtime block rules

- **Always two-step:** literal `object` → `select_first` merge with `runtime_attr_override`.
- Field set is fixed: `cpu_cores`, `mem_gb`, `disk_gb`, `boot_disk_gb`, `preemptible_tries`, `max_retries`, `docker`. Same seven fields, same order, every task.
- `boot_disk_gb: 25` is the default.
- `docker:` images live under `us.gcr.io/broad-dsp-lrma/lr-<name>:<version>` (LRP convention) or a public image (e.g. `staphb/fastqc:0.12.1`). Pin a tag — never `:latest`.
- Disk type: `HDD` by default. `SSD` when IO-bound. `LOCAL` only for very heavy IO (note: LOCAL provisions in 375 GB units).
- `memory:` always suffixed `" GiB"`. `disks:` always prefixed `"local-disk "` and suffixed disk type.
- Column-align the `select_first` calls in `runtime {}` — preserve this when reformatting.

---

## 8. `RuntimeAttr` struct (do not redefine)

```wdl
struct RuntimeAttr {
    Float?  mem_gb
    Int?    cpu_cores
    Int?    disk_gb
    Int?    boot_disk_gb
    Int?    preemptible_tries
    Int?    max_retries
    String? docker
}
```

All optional. Lives in `structs/Structs.wdl`. Imported, never redefined.

---

## 9. Workflow skeleton

```wdl
version 1.0

import "../../../tasks/Utility/PBUtils.wdl"  as PB
import "../../../tasks/Utility/Utils.wdl"    as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ExampleWorkflow {

    meta {
        description: "Multi-sentence overview. What goes in, what comes out, key intermediate steps."

        outputs: {
            final_bam:        "Coordinate-sorted, demultiplexed alignment",
            final_bam_index:  "Index for final_bam"
        }
    }

    parameter_meta {
        ccs_bams:         "GCS path to CCS BAM files"
        ref_map_file:     "table indicating reference sequence and auxiliary file locations"
        gcs_out_root_dir: "GCS bucket to store outputs"
    }

    input {
        Array[File] ccs_bams
        File ref_map_file
        String participant_name
        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ExampleWorkflow/~{participant_name}"

    # ... calls, scatters, conditionals (see §17 for required call-alias convention) ...

    call FF.FinalizeToFile as t_NN_FinalizeBam { input: outdir = outdir, file = t_MM_SomeTask.bam }

    output {
        # re-export key final files / GCS paths; document each in meta.outputs
    }
}
```

### Section order inside `workflow`

1. `meta { ... }` (description; required `outputs:` when workflow has an `output {}` block — see §16)
2. `parameter_meta { ... }` — entry per workflow input
3. `input { ... }`
4. Private declarations (e.g. `Map[String, String] ref_map = read_map(...)`, `String outdir = ...`)
5. `call` / `scatter` / `if` blocks (every call aliased `t_NN_` — see §17)
6. `output { ... }`

### Workflow conventions

- Top-level `outdir`: strip trailing slash with `sub(gcs_out_root_dir, "/$", "")` then append workflow name + participant.
- Call site formatting (call-alias rules in §17):
  - Short call → one line: `call X.Y as t_NN_Name { input: a = b, c = d }`
  - Long call → block form, `input:` on its own indented line, args column-aligned:
    ```wdl
    call PB.Align as t_NN_AlignTranscripts {
        input:
            bam          = t_MM_ClusterTranscripts.clustered_bam,
            ref_fasta    = ref_map['fasta'],
            sample_name  = participant_name,
            map_preset   = "ISOSEQ",
            prefix       = "~{participant_name}.~{BC}",
            runtime_attr_override = { "cpu_cores": 32 }
    }
    ```
- Per-call `runtime_attr_override` is an inline `object` literal — preferred over forking a task to bump CPU.
- Finalization tasks (`FF.FinalizeToFile`, `FF.FinalizeToDir`) live at the **end** of each scatter/branch, with `outdir` paths assembled from the top-level `outdir`.
- Use `select_first([...])` to merge optional outputs from conditional branches.

---

## 10. Naming

| Thing | Convention |
|---|---|
| Task names | `CamelCase` verbs/nouns: `AlignReads`, `ChangeReadGroup`, `MakeChrIntervalList` |
| Workflow names | `CamelCase`, often `<Tech><Description>`: `PBCCSIsoSeq` |
| WDL input vars | `snake_case`: `input_bam`, `ref_fasta`, `gcs_out_root_dir` |
| Output vars | `snake_case`: `aligned_bam`, `aligned_bai` |
| Call aliases (in workflows) | **`t_NN_<CamelCase>`** — see §17 |
| File prefixes | from a `prefix` input, default usually `"out"` |

Read groups / sample fields stay uppercase (`RG`, `ID`, `SM`, `PL`, `LB`) — biology convention, not WDL convention.

---

## 11. Output block

- One declaration per produced file, snake_case.
- Use `glob("...")` for `Array[File]` outputs whose count is dynamic.
- Use `read_map("map.txt")` / `read_tsv(...)` / `read_lines(...)` to surface scalars or tables as typed WDL outputs.
- Every output must be documented in `meta.outputs` (§16).

---

## 12. What to NOT do

- No `version development` or `version 1.1` — repo is on `1.0`.
- No `command { ... }` braces — always heredoc.
- No bare `runtime { docker: "..." cpu: 4 ... }` — always go through `RuntimeAttr`/`select_first`.
- No hardcoded image tags like `:latest`.
- No `bash` outside the heredoc.
- No imports without aliases in pipelines.
- No undocumented inputs (`parameter_meta` must list every input).
- No undocumented outputs (`meta.outputs` must list every output).
- No mixing `${}` and `~{}`. `~{}` is WDL; `${}` is bash. Conflating them silently corrupts substitution.
- No un-aliased `call` statements in workflows — every call uses the `t_NN_` prefix (§17).
- No threading WDL-side `runtime_attr.cpu_cores` / `mem_gb` into the `command` body — derive from `/proc/cpuinfo` and `free -g` (§14).

---

## 13. When reformatting an existing WDL

1. Fix `version 1.0` if missing.
2. Reorder sections per §3 (meta → parameter_meta → input → private decls → command → output → separator → default_attr → runtime_attr → runtime).
3. Convert `command { ... }` to `command <<< ... >>>` and `${}` interpolations to `~{}`.
4. Replace any ad-hoc `runtime { cpu: X memory: Y ... }` with the `RuntimeAttr` two-step + `select_first` pattern (§7).
5. Add `RuntimeAttr? runtime_attr_override` to `input {}` and document in `parameter_meta`.
6. Compute `disk_size` from inputs (§6, default 5×), drop hardcoded `disk_gb` numbers unless tiny static task or override case.
7. Pin docker tags. Move custom-built images to `us.gcr.io/broad-dsp-lrma/lr-<name>:<ver>` only if that's the project's registry — otherwise keep public tag.
8. Add the `#########################` separator before the runtime block.
9. Column-align the runtime `select_first` calls.
10. Insert the §14 resource-detection preamble immediately after `set -euxo pipefail`. Replace hard-coded thread counts and memory flags with `${NUM_CPUS}`, `${JAVA_MEM_GB}`, `${USABLE_RAM_GB}`, `${MEM_PER_THREAD_GB}`. Delete legacy `np=$(...)`, `MEM_FOR_SORT=$(...)` snippets.
11. Add `String extra_args = ""` to `input` block; add corresponding `parameter_meta` entry; append `~{extra_args}` to the primary tool invocation per §15.
12. Add `outputs: { ... }` to `meta {}` block — one entry per `output {}` declaration. Add `tool` / `tool_version` / `tool_url` / `tool_citation` keys to `meta {}` when the task wraps a named tool (§3, §16).
13. **Workflows:** count total `call` statements, pick alias width (2 or 3 digits), rewrite each call's alias to `t_<NN>_<DescriptiveName>` and propagate the rename through all downstream references (§17). Run `miniwdl check <file>` to catch dangling references.

---

## 14. Memory and CPU accounting (mandatory command-block preamble)

**Every `command <<<` block starts with a standard preamble that detects CPUs and RAM, reserves 1 GB for OS, and exports the derived variables.** Compute these even when the current task body doesn't use all of them — keeps tasks reformat-safe.

### Required preamble

Insert immediately after `set -euxo pipefail`:

```bash
set -euxo pipefail

# ---- Resource detection (required preamble) ----
NUM_CPUS=$(grep '^processor' /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')
RAM_IN_GB=$(free -g | grep "^Mem" | awk '{print $2}')

# Reserve 1 GB for OS + container overhead.
USABLE_RAM_GB=$((RAM_IN_GB - 1))
[[ "${USABLE_RAM_GB}" -lt 1 ]] && USABLE_RAM_GB=1

# Per-thread RAM (used by samtools sort -m, etc.)
MEM_PER_THREAD_GB=$(( USABLE_RAM_GB / NUM_CPUS ))
[[ "${MEM_PER_THREAD_GB}" -lt 1 ]] && MEM_PER_THREAD_GB=1

# Java heap (used by GATK/Picard). Same as USABLE_RAM_GB; alias for clarity.
JAVA_MEM_GB=${USABLE_RAM_GB}

echo "NUM_CPUS=${NUM_CPUS}  RAM_IN_GB=${RAM_IN_GB}  USABLE_RAM_GB=${USABLE_RAM_GB}  MEM_PER_THREAD_GB=${MEM_PER_THREAD_GB}  JAVA_MEM_GB=${JAVA_MEM_GB}"
# ---- end preamble ----
```

### Variables produced

| Var | Meaning | Typical consumer |
|---|---|---|
| `NUM_CPUS` | Logical CPUs visible to container | `-@`, `-t`, `--threads`, `-j` |
| `RAM_IN_GB` | Total RAM in container (`free -g`) | Debug/logging |
| `USABLE_RAM_GB` | `RAM_IN_GB - 1`, floored at 1 | Single-process tools |
| `MEM_PER_THREAD_GB` | `USABLE_RAM_GB / NUM_CPUS`, floored at 1 | `samtools sort -m`, `bwa -K` |
| `JAVA_MEM_GB` | Alias of `USABLE_RAM_GB` | `-Xmx`, `-Xms` |

### Usage patterns

**Java tools (GATK, Picard):**
```bash
gatk --java-options "-Xmx${JAVA_MEM_GB}g -Xms${JAVA_MEM_GB}g" HaplotypeCaller \
    -I ~{input_bam} -R ~{ref_fasta} -O ~{prefix}.vcf.gz
```
Set `-Xmx` and `-Xms` equal — avoids GC heap-resize churn.

**samtools sort (per-thread `-m`):**
```bash
samtools sort -@ ${NUM_CPUS} -m ${MEM_PER_THREAD_GB}G ~{input_bam} -o ~{prefix}.bam
```

**Generic multi-threaded tool:**
```bash
some_tool --threads ${NUM_CPUS} --memory ${USABLE_RAM_GB} ...
```

### Rules

- Preamble mandatory in every `command` block — even tiny tasks.
- Never pass full `runtime.memory` to tool. Always use `USABLE_RAM_GB` / `JAVA_MEM_GB`.
- Compute at runtime from `/proc/cpuinfo` and `free -g`. Do **not** thread WDL-side `runtime_attr.cpu_cores` / `mem_gb` into command — stale after `runtime_attr_override`.
- OS reserve: **always 1 GB**. No exceptions.
- Per-thread divides happen **after** OS reserve, never before.
- All floors at 1 GB (avoids `-m 0G`, `-Xmx0g`).

---

## 15. Extra-args passthrough (mandatory input)

**Every task takes an optional `String` input that's appended verbatim to the primary tool invocation.**

### Input declaration

Add to every task's `input` block, immediately before `RuntimeAttr? runtime_attr_override`:

```wdl
input {
    File input_bam
    String prefix

    String extra_args = ""

    RuntimeAttr? runtime_attr_override
}
```

- Name: **`extra_args`**. Always. No synonyms.
- Type: `String` with default `""`. Not `String?`.
- One per task. If a task wraps multiple distinct tools (rare — split it), only the **primary** tool gets `extra_args`.

### `parameter_meta` entry

```wdl
parameter_meta {
    input_bam:             "..."
    prefix:                "..."
    extra_args:            "Additional command-line args appended verbatim to the samtools sort invocation"
    runtime_attr_override: "Override the default runtime attributes"
}
```

Replace `samtools sort` with the actual binary name.

### Usage in `command`

Append `~{extra_args}` to the **primary** tool invocation, after the structured flags:

```bash
samtools sort -@ ${NUM_CPUS} -m ${MEM_PER_THREAD_GB}G ~{extra_args} ~{input_bam} -o ~{prefix}.bam
```

```bash
gatk --java-options "-Xmx${JAVA_MEM_GB}g -Xms${JAVA_MEM_GB}g" HaplotypeCaller \
    -I ~{input_bam} \
    -R ~{ref_fasta} \
    -O ~{prefix}.vcf.gz \
    ~{extra_args}
```

```bash
minimap2 -ayYL --MD --eqx -x ~{map_preset} -t ${NUM_CPUS} ~{extra_args} ~{ref_fasta} ~{reads} > tmp.sam
```

### Rules

- `~{extra_args}` placement: **before positional arguments**, **after** the structured flags the task controls. Lets users override defaults via tool's last-flag-wins semantics.
- Do **not** quote `~{extra_args}` in bash — that collapses multi-flag strings into one argv element. Plain `~{extra_args}` lets the shell tokenize.
- Do **not** validate or sanitize. Pass through verbatim.
- If task has multiple sequential tool calls (e.g. sort + index), `extra_args` applies to the primary invocation only. Document in `parameter_meta`.
- Default `""` is safe — empty interpolation drops to nothing in bash.

---

## 16. Output documentation

WDL 1.0 has no formal output `parameter_meta`. Use the `meta` block with a free-form `outputs` key. Engines (Cromwell, miniwdl) preserve arbitrary `meta` keys.

### Required: `meta.outputs`

Every task and workflow with an `output {}` block must document each output:

```wdl
meta {
    description: "..."
    tool: "samtools sort"
    tool_version: "1.23.1"
    tool_url: "https://www.htslib.org/"

    outputs: {
        output_bam: "Coordinate-sorted BAM",
        output_bai: "BAM index (.bam.bai)"
    }
}
```

### Rules

- Key: **`outputs`**. Always. Inside `meta {}`, not a separate top-level block.
- Value: object literal with one key per `output {}` declaration. Key names must **match** the `output` variable names exactly.
- Value type: short string description. Object form allowed for tagged outputs but strings preferred.
- Order: same order as the `output {}` block.
- If output is a glob result, document the glob pattern in the description: `"All per-contig interval files matching contig.*.intervals"`.
- Workflows: document **only outputs declared in the workflow's own `output {}` block.** Don't re-document task outputs.

### Engine note

Cromwell and miniwdl ignore unknown `meta` keys at execution. They surface them in `womtool inputs` / `miniwdl check` JSON output, so downstream doc generators can scrape them.

---

## 17. Workflow call aliasing (mandatory)

**Every `call` in a workflow file is aliased.** Alias prefix is `t_NN_` (or `t_NNN_`), where the number is the call's position in lexical source order.

### Format

```
t_<NN>_<DescriptiveName>
```

- `<NN>` = zero-padded ordinal of the `call` statement in the file, counting from top to bottom in source order, starting at `01`.
- Width:
  - `2` digits (`t_01_`, `t_99_`) when the workflow contains **fewer than 100** `call` statements.
  - `3` digits (`t_001_`, `t_999_`) when ≥ 100.
  - Pick width once per file based on total count; do not mix widths in one workflow.
- `<DescriptiveName>` = CamelCase, describes purpose at this call site. Mirrors §10 alias rules.

### Examples

```wdl
call Utils.MergeBams as t_01_MergeAllReads { input: bams = ccs_bams, prefix = sample }

call PB.Demultiplex as t_02_Demultiplex {
    input:
        bam          = t_01_MergeAllReads.merged_bam,
        prefix       = sample,
        barcode_file = barcode_file
}

scatter (demux_bam in t_02_Demultiplex.demux_bams) {
    call PB.RefineTranscriptReads as t_03_RefineTranscriptReads { input: bam = demux_bam }
    call PB.ClusterTranscripts    as t_04_ClusterTranscripts    { input: bam = t_03_RefineTranscriptReads.refined_bam }
    call PB.Align                 as t_05_AlignTranscripts      { input: bam = t_04_ClusterTranscripts.clustered_bam }
}

call FF.FinalizeToFile as t_06_FinalizeBam { input: outdir = outdir, file = t_05_AlignTranscripts.aligned_bam }
```

### Rules

- **Every `call` gets an alias.** No exception, even single-use calls and calls to differently-named tasks.
- **Numbering is by `call` statement, not runtime invocation.** A single `call` inside a `scatter` still consumes one number, regardless of how many shards execute.
- **Order is source-file lexical order, top to bottom.** Calls nested inside `scatter`, `if`, or sub-block count in the order they textually appear.
- **Sub-workflow calls (`call SubWorkflow as ...`)** count too — same `t_NN_` prefix.
- **Don't renumber on minor edits.** When inserting a new call mid-file, renumber every subsequent call. This is intentional friction; the renumbering keeps order monotonic and `t_NN_` references in `input:` blocks consistent.
- **Downstream references update with the alias.** Renaming `t_03_Foo` → `t_04_Foo` requires updating every `t_03_Foo.<output>` reference. Use editor find-replace; verify with `miniwdl check`.
- **Workflow outputs that re-export task outputs** also use the `t_NN_` alias: `File final_bam = t_05_AlignTranscripts.aligned_bam`.
