# Project Tree

```
.
├── API.md
├── bridge
│   ├── ancestry_simulation.py
│   ├── configuration.py
│   ├── demography_builder.py
│   ├── header_dataclasses.py
│   ├── __init__.py
│   ├── loci_parser.py
│   ├── observed_data.py
│   ├── parameter_sampling.py
│   ├── pipeline.py
│   ├── prior_parser.py
│   ├── __pycache__
│   ├── reftable_loop.py
│   ├── scenario_parser.py
│   ├── snp_writer.py
│   ├── statobs_parser.py
│   ├── stats_group_parser.py
│   └── summary_statistics.py
├── CLAUDE.md
├── docs
├── notebook
│   ├── comparaison_summary.csv
│   ├── comparaison_summary_exemple5.csv
│   ├── comparaison_summary_fixed_params.csv
│   ├── comparaison_summary_human.csv
│   ├── compare_raw_genotypes_fixed_params.ipynb
│   ├── compare_reftables_exemple2.ipynb
│   ├── compare_reftables_exemple5.ipynb
│   ├── compare_reftables_human.ipynb
│   └── test_fonction.ipynb
├── notes
│   ├── api.md
│   ├── commits.md
│   ├── exploration.md
│   ├── report.md
│   ├── resume_coalescence_lignees.md
│   ├── resume_stat_dna.md
│   └── tree.md
├── pyproject.toml
├── README.md
├── reference
│   ├── human
│   ├── toy_example1_ms
│   ├── toy_example2_ms_dna
│   ├── toy_example2_ms_dna_50loci
│   ├── toy_example2_ms_dna_XY
│   ├── toy_example3
│   ├── toy_example3_500loci
│   ├── toy_example4
│   ├── toy_example4_MRC1
│   ├── toy_example5
│   └── toy_example5_500loci
├── scripts
│   ├── benchmark_1000.py
│   ├── calibrate_reftable.py
│   ├── compare_reftable_te4.py
│   ├── generate_test_reftable.py
│   ├── param_keepers.py
│   ├── priors_keeper.py
│   ├── profile_one_particle.py
│   ├── profile_python_stats.py
│   ├── profile_sim.py
│   ├── profile_stats.py
│   ├── __pycache__
│   ├── replay_diyabc_priors_dna_50loci.py
│   ├── replay_diyabc_priors_dna.py
│   ├── replay_diyabc_priors.py
│   ├── run_test.py
│   ├── validate_observed_stats_poolseq.py
│   └── validate_stats.py
├── tests
│   ├── conftest.py
│   ├── __pycache__
│   ├── test_ancestry_simulation.py
│   ├── test_demography_builder.py
│   ├── test_loci_parser.py
│   ├── test_observed_data.py
│   ├── test_parameter_sampling.py
│   ├── test_pipeline.py
│   ├── test_prior_parser.py
│   ├── test_reftable_loop.py
│   ├── test_scenario_parser.py
│   ├── test_snp_writer.py
│   ├── test_stats_group_parser.py
│   └── test_summary_statistics.py
├── tmp
│   ├── benchmark_1000
│   ├── human_scratch_test
│   ├── profile_one
│   ├── replay_diyabc_priors
│   ├── validate_batching
│   └── validate_stats
└── tools
    ├── generate_api_md.py
    └── generate_report.py

29 directories, 67 files

```