# calculate_psi

Count splice-junction reads supporting each splicing isoform from STAR `SJ.out.tab` files, and calculate the ratio (PSI) of each isoform.  
Any number of isoforms (N ≥ 2) per gene is supported.

## Requirements

- Python 3
- [pysam](https://github.com/pysam-developers/pysam)
- [HTSlib](https://github.com/samtools/htslib) (`bgzip`, `tabix`), for preparing input files

## Isoform examples

|====| : exon, ----- : junction (intron spliced out)

**FGFR2 IIIb / IIIc (mutually exclusive exons, 2 junctions per isoform)**

```
           A              IIIc           IIIb              B
2c    |=========|-------|======|----------------------|=========|
2b    |=========|----------------------|======|-------|=========|
```


**GENE_X (3 isoforms, 1 junction per isoform)**

```
        E1              E2          E3          E4
iso1  |====|----------|====|
iso2  |====|----------------------|====|
iso3  |====|----------------------------------|====|
```

## Input

### 1. SJ.out.tab file (bgzip-compressed and tabix-indexed)

Junction files produced by STAR (`SJ.out.tab`). Each file must be sorted, bgzip-compressed and tabix-indexed.

```bash
sort -k1,1 -k2,2n SJ.out.tab | bgzip > SAMPLE.SJ.out.tab.gz
tabix -s 1 -b 2 -e 3 SAMPLE.SJ.out.tab.gz    # creates SAMPLE.SJ.out.tab.gz.tbi
```


### 2. Sample list (`-input_list`)

Tab-separated, no header.

| Column | Description |
|---|---|
| 1 | Sample name |
| 2 | Path to the bgzip-compressed, tabix-indexed SJ.out.tab file |

```
SAMPLE_A	/path/to/SAMPLE_A.SJ.out.tab.gz
SAMPLE_B	/path/to/SAMPLE_B.SJ.out.tab.gz
```

### 3. Isoform junction list (`-sj_file`)

Tab-separated, no header.  
One gene per line.  
The number of isoforms per gene ≥ 2

| Gene | junctions of isoform 1 | … | junctions of isoform N | name of isoform 1 | … | name of isoform N |
|---|---|---|---|---|---|---|

- Junctions are written as `chr:start-end` (same coordinate system as `SJ.out.tab`: first and last intronic base, 1-based).
- If an isoform is defined by multiple junctions, separate them with commas (`,`). Reads of all junctions of an isoform are summed.

```
# Example #
GENE_X	chr1:1000-2000	chr1:1000-2500	chr1:1000-3000	iso1	iso2	iso3
```
```
# Example #
FGFR2	chr10:121515320-121517318,chr10:121517464-121519978	chr10:121515320-121518681,chr10:121518830-121519978	2c	2b
```

## Usage

```bash
python calculate_psi.py \
    -input_list sample_list.txt \
    -sj_file sj_list.txt \
    -output_prefix result/psi_result_
```

| Option | Description | Default |
|---|---|---|
| `-input_file` | Sample list (see above) | (required) |
| `-sj_file` | Isoform junction list (see above) | (required) |
| `-output_prefix` | Prefix of output files | `psi_result_` |

## Output

One tab-separated file per gene: `<output_prefix><Gene>.txt`

| Column | Description |
|---|---|
| `Gene` | Gene name |
| `Sample` | Sample name |
| `Isoforms` | Isoform names joined by `\|` |
| `<isoform>_count` | Read count of each isoform (N columns) |
| `SUM(isoforms)` | Total read count of all isoforms |
| `<isoform>_ratio` | Count of each isoform / `SUM(isoforms)` (N columns) |

Example (`psi_result_FGFR2.txt`):

```
Gene	Sample	Isoforms	2c_count	2b_count	SUM(isoforms)	2c_ratio	2b_ratio
FGFR2	SAMPLE_A	2c|2b	0	11	11	0.0	1.0
FGFR2	SAMPLE_B	2c|2b	136	114	250	0.544	0.456
FGFR2	SAMPLE_C	2c|2b	33	266	299	0.110	0.890
```

### Notes on ratios

- Ratios are rounded to 3 decimal places.
- `0.0` and `1.0` are reported only when the ratio is exactly 0 or 1. Values that would round to 0 or 1 are reported as `0.001` / `0.999`.
- If `SUM(isoforms)` is 0, all ratios are reported as `0.0`.
- Ratios from samples with a small `SUM(isoforms)` are less reliable; check the count columns when interpreting them.

## Author

Naoko Iida

## License

MIT

## Cite
