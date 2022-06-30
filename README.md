### Overview

### How to run

The CHTC serpent workflow manager, developed by David A. Baker, allows extremely precise, case-specific configuration with JSON files. In this repo, these config files are located in `iwes_genotyping/config`.

_NOTE:_ Only config files that you intend to use should be in `iwes_genotyping/config`. Additional config files will confuse the workflow with conflicting directions.

That said, here's how to genotype iWES FASTQs using this repo:
1. Edit the following JSON-formatted configuration files with _absolute file paths_:
	- ``
	-
	-
2. Change into the project working directory (the directory that you cloned from [the iWES-genotyper GitHub repo](https://github.com/nrminor/iWES-genotyper)).
3. Run the following command:
	```
	python3 chtc_serpent_v2/main.py \
		--config_dir iwes_genotyping/config \
		--submission_dir ~/chtc_serpent/chtc_serpent_submissions
	```

