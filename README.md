# QCxMS
[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fanie.201300158%20-blue)](https://doi.org/10.1002/anie.201300158) [![DOI](https://img.shields.io/badge/DOI-10.1021%2Facsomega.9b02011%20-blue)](https://doi.org/10.1021/acsomega.9b02011)

This is the download repository for the QCxMS program. Currently, only the executeable is provided on the [latest release page](https://github.com/qcxms/QCxMS/releases/tag/latest).

**Installation**

Untar the zipped archive:

```bash
tar -xvzf QCxMS.vX.X.tar.xz
```

The following files are being extracted: `qceims` `pqceims` `q-batch` `getres` `.mass_raw.agr` `.XTBPARAM` `EXAMPLE`

Place the executables into your ``$HOME/bin/`` directory or path. Place the `.XTBPARAM` folder and `.mass_raw.arg` file into your `$HOME` directory (these files can appear to be hidden). 
For  use the [PlotMS](https://github.com/qcxms/PlotMS) program. The of the calculated spectra we recommend the useage of ``xmgrace``. Examples for an input file and coordinate file are to be found in the `EXAMPLES` file.


**Documentation**

For a full documentation on installation and input settings can be fond at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/qcxms_doc/qcxms.html).



**From QCEIMS to QCxMS:**
- All names have been changed from `qceims.xxx` to `qcxms.xxx`.
- The `q-batch`, `pqcxms` and `plotms` script have been updated.
- Collision induced dissociation (CID) calculations are now available. Set *cid* in the `qcxms.in` file. 
- GFN1- and GFN2-xTB calculations can be used for CID calculations. DFT methods are not yet tested.

