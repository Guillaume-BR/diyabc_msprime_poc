# Project Tree

```
.
├── benchmark_1000.py
├── bridge
│   ├── ancestry_simulation.py
│   ├── demography_builder.py
│   ├── __init__.py
│   ├── loci_parser.py
│   ├── observed_data.py
│   ├── parameter_sampling.py
│   ├── pipeline.py
│   ├── prior_parser.py
│   ├── __pycache__
│   │   ├── ancestry_simulation.cpython-311.pyc
│   │   ├── demography_builder.cpython-311.pyc
│   │   ├── __init__.cpython-311.pyc
│   │   ├── loci_parser.cpython-311.pyc
│   │   ├── observed_data.cpython-311.pyc
│   │   ├── parameter_sampling.cpython-311.pyc
│   │   ├── pipeline.cpython-311.pyc
│   │   ├── prior_parser.cpython-311.pyc
│   │   ├── reftable_loop.cpython-311.pyc
│   │   ├── scenario_parser.cpython-311.pyc
│   │   ├── scenario_types.cpython-311.pyc
│   │   ├── snp_writer.cpython-311.pyc
│   │   ├── statobs_parser.cpython-311.pyc
│   │   └── summary_statistics.cpython-311.pyc
│   ├── reftable_loop.py
│   ├── scenario_parser.py
│   ├── scenario_types.py
│   ├── snp_writer.py
│   ├── statobs_parser.py
│   └── summary_statistics.py
├── calibrate_reftable.py
├── docs
├── generate_test_reftable.py
├── msprime_cpp
│   ├── msprime_from_cpp
│   └── msprime_from_cpp.cpp
├── notes
│   ├── api.md
│   ├── commits.md
│   ├── exploration.md
│   ├── report.md
│   └── tree.md
├── param_keepers.py
├── priors_keeper.py
├── profile_one_particle.py
├── profile_python_stats.py
├── profile_sim.py
├── profile_stats.py
├── pyproject.toml
├── README.md
├── reference
│   └── human
│       ├── headerRF.txt
│       ├── header.txt
│       ├── human_snp_all22chr_maf5.snp
│       ├── reftableRF.bin
│       └── RNG_state_0000.bin
├── run_test.py
├── tests
│   ├── __pycache__
│   │   └── test_scenario1_human.cpython-311-pytest-9.0.3.pyc
│   └── test_scenario1_human.py
├── tmp
│   ├── bench_reftable.bin
│   └── validate_stats
│       ├── first_records_of_the_reference_table_0.txt
│       ├── headerRF.txt
│       ├── header.txt
│       ├── human_snp_all22chr_maf5.snp
│       ├── human_snp_all22chr_maf5.snp.bin
│       ├── human_snp_all22chr_maf5.snpbin.txt
│       ├── maf.txt
│       ├── reftable.log
│       ├── reftableRF.bin
│       ├── RNG_state_0000.bin
│       └── statobsRF.txt
├── tools
│   ├── generate_api_md.py
│   └── generate_report.py
└── validate_stats.py

12 directories, 69 files

```