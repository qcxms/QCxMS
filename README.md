# QCxMS
[![License](https://img.shields.io/github/license/qcxms/qcxms)](https://github.com/grimme-lab/xtb/blob/main/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/qcxms/qcxms)](https://github.com/qcxms/QCxMS/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fanie.201300158%20-blue)](https://doi.org/10.1002/anie.201300158) [![DOI](https://img.shields.io/badge/DOI-10.1021%2Facsomega.9b02011%20-blue)](https://doi.org/10.1021/acsomega.9b02011)

This is the download repository for the QCxMS program. 

**Installation**

### Binary 

Statically linked binaries (Intel Compiler 21.3.0) can be found at the [latest release page](https://github.com/qcxms/QCxMS/releases/latest).

Untar the zipped archive:

```bash
tar -xvzf QCxMS.vX.X.tar.xz
```

The following files are being extracted: `qcxms` `pqcxms` `q-batch` `getres` `.XTBPARAM` `EXAMPLE`

Place the executables into your ``$HOME/bin/`` directory or path. Place the `.XTBPARAM` folder and `.mass_raw.arg` file into your `$HOME` directory (these files can appear to be hidden). 

### Conda

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/qcxms.svg)](https://anaconda.org/conda-forge/qcxms)

Installing `qcxms` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
```

Once the `conda-forge` channel has been enabled, `qcxms` can be installed with:

```
conda install qcxms
```

It is possible to list all of the versions of `qcxms` available on your platform with:

```
conda search qcxms --channel conda-forge
```


### Meson

Using [meson](https://mesonbuild.com/) as build system requires you to install a fairly new version like 0.57.2 or newer.
To use the default backend of meson you have to install [ninja](https://ninja-build.org/) version 1.10 or newer.

```bash
export FC=ifort CC=icc
meson setup build -Dfortran_link_args=-static
ninja -C build 
```

This will build a static linked binary in the ``build`` folder. Copy the binary from ``build/qcxms`` file into a directory in your path, e.g. ``~/bin/``.


**Documentation**

A more detailed documentation on topics like input settings can be fond at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms.html). 
Examples to test QCxMS can be found in the `EXAMPLES` folder. Here, input and coordinate files are provided for either EI or CID run modes. 


**From QCEIMS to QCxMS:**
- All names have been changed from `qceims.xxx` to `qcxms.xxx`.
- The `q-batch`, `pqcxms` and `plotms` script have been updated.
- Collision induced dissociation (CID) calculations are now available. Set *cid* in the `qcxms.in` file (see
  documentation) 

**The tblite library for xTB calculations**
- The [tblite](https://github.com/awvwgk/tblite) library has been included into the program code. This keeps xtb up-to-date and decreases the computational time for calculations done with GFN1- and GFN2-xTB when compared to earlier versions. 


**Plotting Spectra**

To evaluate the results and create a spectrum, download and use the [PlotMS](https://github.com/qcxms/PlotMS) program. 
The [documentation](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms_plot.html) explains the basic
functionalities of the program. 

The program provides *mass.agr*, *JCAMP-DX* and *.csv* are files that can be analyzed. 
For visualization of the calculated spectra, we recommend the usage of the **xmgrace** program. 

### Updates

Versions PlotMS v.6.0 and higher now provide **exact masses**.
Experimental files in `.csv` format can now be read and plotted against the computed spectra.
The `.mass_raw.agr` file was moved to the [PlotMS](https://github.com/qcxms/PlotMS) repository. 
