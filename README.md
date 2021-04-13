# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.
You can fork this repository to build your own pipeline.
To get the last commits from this repository into your fork use the following commands:

For the first time:
```sh
git remote add upstream git@gitbio.ens-lyon.fr::pipelines/nextflow.git
git pull upstream master
```

Then to make an update:
```sh
git pull upstream master
git merge upstream/master
```

**If you created your `.config` file before version `0.4.0` you need to run the script `src/.update_config.sh` to use the latest docker, singularity and conda configuration (don't forget to check your config files afterward for typos).**

## Getting Started

These instructions will get you a copy of the project as a template when you want to build your own pipeline.

[you can follow them here.](doc/getting_started.md)

## Contributing

If you want to add more tools to this project, please read the [CONTRIBUTING.md](CONTRIBUTING.md).

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitbio.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details
