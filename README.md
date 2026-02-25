[![Rust](https://github.com/codialab/sqz/actions/workflows/rust.yml/badge.svg)](https://github.com/codialab/sqz/actions/workflows/rust.yml)

# Human-readable path compression in GFA-formatted pangenome graphs

<p align="center">
  <img src="logo.png" height=500px />
</p>


`sqz` is a tool for compressing and decompressing path annotations in GFA files. It uses a modified version of the GFA-format containing `Q`- and `Z`-lines.

## Usage

Use
```bash
sqz compress-full <GFA-FILE> > <COMPRESSED-GFA-FILE>
```
to compress a GFA file and
```bash
sqz decompress <COMPRESSED-GFA-FILE> > <GFA-FILE>
```
for decompression.
For large GFA files take a look at
```bash
sqz compress-partial --help
```

## Installation
### Using a precompiled binary
Use
```bash
wget 
chmod +x 
```
and make `sqz` available on your `PATH` to install it.

### Building from source
`sqz` is written in [RUST](https://rust-lang.org) and requires a working RUST build system (version >= 1.74.1) for installation. See [here](https://www.rust-lang.org/tools/install) for more details.

```bash
git clone git@github.com:codialab/sqz.git
cd sqz
cargo build --release
```

## Format
The format of `sqz` is mostly based on the [GFA format](https://gfa-spec.github.io/GFA-spec/GFA1.html), but extended to include one new types of lines: `Q`-lines.
`W`-lines have been changed to allow the usage of identifiers from `Q`-lines.

### `Q` Rule line

A Q-line defines part of a compressed walk that can be used as part of other
compressed walks.

#### Required fields

| Column | Field             | Type      | Regexp              | Description
|--------|-------------------|-----------|---------------------|------------
| 1      | `RecordType`      | Character | `Q`                 | Record type
| 2      | `Name`            | String    | `@[!-)+-<>-~][!-~]*` | Rule name
| 3      | `CompressedWalk`  | String    | `([><][!-;=?-~]+)+` | Compressed Walk

A `Walk` is defined as
```txt
<walk> ::= ( `>' | `<' <segId> )+
```
where `<segId>` corresponds either to the identifier of a segment or the
identifier of a Q-line. A valid walk must exist in the graph. The identifier of
a Q-line starts with the character `@`


### `W` Compressed walk line

A walk line describes an oriented walk in the graph. It is only intended for a
graph without overlaps between segments.
Note that W-lines can not use jump connections (introduced in v1.2).

#### Required fields

| Column | Field             | Type      | Regexp                   | Description
|--------|-------------------|-----------|--------------------------|------------
| 1      | `RecordType`      | Character | `W`                      | Record type
| 2      | `SampleId`        | String    | `[!-)+-<>-~][!-~]*`      | Sample identifier
| 3      | `HapIndex`        | Integer   | `[0-9]+`                 | Haplotype index
| 4      | `SeqId`           | String    | `[!-)+-<>-~][!-~]*`      | Sequence identifier
| 5      | `SeqStart`        | Integer   | `\*\|[0-9]+`             | Optional Start position
| 6      | `SeqEnd`          | Integer   | `\*\|[0-9]+`             | Optional End position (BED-like half-close-half-open)
| 7      | `CompressedWalk`  | String    | `([><][!-;=?-~]+)+`      | Compressed Walk

For a haploid sample, `HapIndex` takes 0. For a diploid or polyploid sample,
`HapIndex` starts with 1. For two W-lines with the same
(`SampleId`,`HapIndex`,`SeqId`), their [`SeqSart`,`SeqEnd`) should have no
overlaps. A `Walk` is defined as
```txt
<walk> ::= ( `>' | `<' <segId> )+
```
where `<segId>` corresponds either to the identifier of a segment or the
identifier oft a Q-line. A valid walk must exist in the graph.

### Example

```txt
S       s11     ACCTT
S       s12     TC
S       s13     GATT
L       s11     +       s12     -       0M
L       s12     -       s13     +       0M
L       s11     +       s13     +       0M
Q       @q1      >s11<s12
W       NA12878 1       chr1    0       11      >@q1>s13
