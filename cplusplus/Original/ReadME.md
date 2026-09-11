# `cplusplus/Original/` Folder Structure

The multithreaded framework implemented here is based on the island PSO model described in:

```bibtex
@article{abadlia2022island,
  author  = {Abadlia, Houda and Smairi, Nadia and Ghedira, Khaled},
  title   = {Comparative performance evaluation of island particle swarm algorithm applied to solve constrained and unconstrained optimization problems},
  journal = {Journal of Intelligent \& Fuzzy Systems},
  volume  = {43},
  pages   = {1-17},
  year    = {2022},
  month   = {03},
  doi     = {10.3233/JIFS-213380}
}
```

| Subfolder | Description |
|---|---|
| `code_ideas/` | Preliminary/experimental implementations of the multithreaded framework and the l-back functionality |
| `converters/` | Utility Python scripts (`plu_to_truth_table.py`, `truth_table_to_plu.py`, `recursive_parser.py`) for converting between truth table and packed PLU formats. These will later be integrated into a web application pipeline as part of the research methodology. Also contains `plu_files/`, the PLU-format truth table files |
| `csvs/` | CSV output files logging evolution results from experiment runs |
| `data_files/` | `.dta` input files with truth tables in packed PLU format, organized into `demux/`, `full_adder/`, `icomp/`, and `tests/` (the original Coello `.dta` files, converted to packed PLU representation) |
| `results/` | **Note:** contains results from a *previous* iteration of the code — prior to the multithreaded rewrite and before l-back functionality was added. Kept for reference/comparison only |

## Compilation

### `cplusplus/Original/` (Current Implementation)

This folder includes a Makefile. To build:

```bash
make
```

To clean the build:

```bash
make clean
```