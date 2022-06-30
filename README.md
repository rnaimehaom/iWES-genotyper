## Overview

## Installation

To install the workflow along with its submodules, simply run the following command in your directory of choice:
	```
	git clone --recursive https://github.com/nrminor/iWES-genotyper.git .
	```

The `--recursive` flag is necessary to pull the git submodule [serpent](https://github.com/DABAKER165/chtc_serpent_v2) along with the rest of the repo.

## How to run

The CHTC serpent workflow manager, developed by [David A. Baker](https://github.com/DABAKER165), allows extremely precise, case-specific configuration with JSON files. In this repo, these config files are located in `iwes_genotyping/config`.

_NOTE:_ Only config files that you intend to use should be in `iwes_genotyping/config`. Additional config files will confuse the workflow with conflicting directions.

That said, here's how to genotype iWES FASTQs using this repo:
1. Edit the following JSON-formatted configuration files with local and remote filepaths:
	- `iwes_genotyping/config/default_iwes_artic.json` - This file configures how the genotyping workflows run and where they run. File paths for your local machine and for your CHTC directory must be specified. The settings also specify a Google Drive output directory and an data backup array local to groups at the [AIDS Vaccine Research Laboratory](https://dholk.primate.wisc.edu/wiki/home/page.view?name=home_index). The file paths in this config file should be _absolute_.
	- `iwes_genotyping/config/chtc_iwes_settings.json` - If you are running this workflow from within the [O'Connor Group](https://github.com/dholab), we do not recommend you change these settings. Otherwise, you will need to change the values for "priority_flag", or remove this setting entirely.
	- `iwes_genotyping/config/config.yaml` - This YAML-formatted config file is only used by the snakemake workflows in `iwes_genotyping/static_files`. As such, the file paths in this config file should be _relative_. If you re-organize this repo, you will not need to change this config file.
2. Change into the project working directory (the directory that you cloned from [the iWES-genotyper GitHub repo](https://github.com/nrminor/iWES-genotyper)).
3. Run the following command:
	```
	python3 chtc_serpent_v2/main.py \
		--config_dir iwes_genotyping/config \
		--submission_dir ~/chtc_serpent/chtc_serpent_submissions
	```

