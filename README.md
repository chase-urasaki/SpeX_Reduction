# Spectral Reduction Pipeline

## Overview
This repository contains code for reducing astronomical spectra. The pipeline will handle essential steps in spectral reduction, including:

- Bias and dark subtraction
- Flat-field correction
- Wavelength calibration
- Flux calibration
- Sky subtraction
- Combination of spectral orders (if applicable)

## Features
- Support on IRTF/SpeX

## Requirements
- Python 3.8+
- Dependencies:
  - numpy
  - scipy
  - matplotlib
  - astropy
  - specutils (optional for advanced spectral analysis)
  - pandas (optional for data management)

## Installation
Clone this repository and install the required dependencies:

```bash
git clone https://github.com/chase-urasaki/SpeX_Reduction.git
pip install -r requirements.txt
```

## Usage

### Basic Example
To process a set of raw spectral files:

```bash
python reduce_spectra.py --input /path/to/raw/files --config config.yaml --output /path/to/output
```

### Configuration
The pipeline uses a YAML configuration file to specify parameters such as:
- Input data format
- Calibration file locations
- Reduction steps to apply

See `config_example.yaml` for a template.

## File Structure
```
repository/
├── data/                  # Example data files
├── docs/                  # Documentation
├── scripts/               # Core reduction scripts
├── tests/                 # Unit tests
├── config_example.yaml    # Example configuration file
├── requirements.txt       # Python dependencies
└── README.md              # This file
```

## Contributing
Contributions are welcome! Please:
1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Submit a pull request with a clear description of your changes

## License

This project was developed as part of ASTR633 - Astrophysical Techniques
