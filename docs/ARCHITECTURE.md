# Architecture

## Package Layout

- `src/methphaser/cli.py`: minimal CLI (`methphaser`)
- `src/methphaser/pipeline.py`: one-shot orchestration entry
- `src/methphaser/discovery.py`: GTF/VCF auto-discovery from BAM context
- `src/methphaser/step1_core.py`: methylation assignment and block-relationship inference core
- `src/methphaser/step1_parallel.py`: full-GTF chromosome/block parallel execution
- `src/methphaser/step1_regions.py`: phase-region planning (legacy full-span and boundary-window modes)
- `src/methphaser/step2_postprocess.py`: altered VCF/BAM generation from relationships

## Command Surface

Only one command is exposed:

- `methphaser`

## Design Decisions

- Truth-VCF benchmarking code paths removed from runtime flow.
- Full-span phase regions are the default for legacy parity.
- Boundary-window phase regions are available as a faster optional mode.
- BAM merge is mandatory and handled automatically in the pipeline.

## Tests

- `tests/test_discovery.py`: discovery logic
- `tests/test_step1_parallel_helpers.py`: step1 helper behavior
- `tests/test_integration_smoke.py`: end-to-end smoke test on local subset data in `tests/data`
