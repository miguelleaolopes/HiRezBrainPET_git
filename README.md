# CASToR Reconstruction for HiRezBrainPET

[CASToR](https://castor-project.org/) is an open source software for tomographic reconstruction. This repository contains modifications to the code and new the scripts to run the reconstruction for the HiRezBrainPET system, an innovative PET scanner that uses RPC as detectors. The [first results](https://arxiv.org/abs/2211.05860) were published in November 2022 by Paulo Fonte et al. using a simple reconstruction process, and the goal with this repository is to make it more generic and customizable to the system, and allowing a more flexible and scalable tool for researchers in the field of medical imaging, or even the broader public, to use and replicate the results.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## License

[MIT](https://choosealicense.com/licenses/mit/)
