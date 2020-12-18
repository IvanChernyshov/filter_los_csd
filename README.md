# filter_los_csd

This is the manual for the script which filters line-of-sight contacts from the results of ConQuest search over the Cambridge Structure Database (CSD) entries. The script is available via [GitHub](https://github.com/IvanChernyshov/filter_los_csd). For the underlying scientific context see the corresponding work: DOI [10.1002/cphc.201901083](http://doi.org/10.1002/cphc.201901083).

If you use the script in your work, please cite generously!

I. Yu. Chernyshov, I. V. Ananyev, E.A. Pidko, Revisiting van der Waals Radii: From Comprehensive Structural Analysis to Knowledge‐Based Classification of Interatomic Contacts. _ChemPhysChem_ **2020**, _21_, 370–376, DOI: [10.1002/cphc.201901083](http://doi.org/10.1002/cphc.201901083).

## Intro

This manual consists of two parts, describing installation details and syntax of the script. The example of input and output files for the script is provided in [*example*](https://github.com/EPiCs-group/filter_los_csd/tree/master/example) directory.

## Installation

*filter_los_csd* is a multiplatform Python3 package. Use `pip` to install it:

```
> pip install filter_los_csd
```

*filter_los_csd* depends on the [PyCifRW](http://pypi.org/project/PyCifRW/) package, which requires a C/C++ compiler, and pip installation will fail if you do not have one. The best choice for Windows is `Visual Studio Building Tools`, and for Linux, `gcc` would be enough. If you forget about it, pip will give you the installation error and platform-specific advice on fixing it.

A `Permission denied` error may occur when running the script in some versions of Linux:

```
> filter_los_csd.py input.csv input.cif
-bash: /path/to/script/filter_los_csd.py: Permission denied
```

To fix it, set `umask` to `022` before script installation (see [ref](https://stackoverflow.com/questions/36898474/how-to-install-a-module-for-all-users-with-pip-on-linux) for more details), or just change permissions:

```
> sudo chmod 755 /path/to/script/filter_los_csd.py
```

## Syntax

### Input and Output

The script takes two necessary parameters as input:

* *path_csv*: path to CSV file, containing info on parameters of contacts in question.
  **Please note, that CSV file must contain labels of atoms forming the contact and the contact distance!**
* *path_cif*: path to multiple CIF file, containing CSD entries from CSV file. **Please note, that CIF file must contain info on bonds.** Make sure that the corresponding checkbox was selected before downloading CIF.
  **Also note, that script were tested only for CIF files extracted from CSD (ConQuest)!**

Thus, the easiest command is:

```
> filter_csd_los.py test.csv test.cif
```

As output the script creates **{csv_name}_los.csv** file containing the same info as the original file with three additional columns:

* "LOS": has four possible values:
  * "+": the corresponding contact is line-of-sight;
  * "–": the corresponding contact is not line-of-sight;
  * "?": the corresponding crystal is larger than specified cell volume cutoff;
  * "!": the error occurred during the calculation. The only type of errors caught during the testing is inability to find a contact with given atomic labels and distance. The main source of these errors are disorder issues, thus tuning *--tol* parameter can solve the problem.
* "SHIELDING": contact shielding value (the definition of the term is given in the manuscript);
* "SHIELD_ATOM": label of the shielding atom.

The contact and the shielding atom can be visualized in two steps:

1. find the contact by atomic labels and distance;
2. the nearest atom to the contact line with SHIELD_ATOM label is a shielding atom.

### Optional parameters

In addition, the script has several optional parameters:

* *-r* or *--radii*: type of van der Waals radii used for shielding calculation, "csd" by default. Available values are:
  * "csd": version usually used in CCDC products (ConQuest, Mercury, etc). It is same as Bondi's version but with r(H) = 1.09 Å; unknown atomic radii was set to 2.0 Å;
  * "bondi": Bondi version; unknown atomic radii set to 2.0 Å;
  * "rt": Rowland&Taylor version, unknown atomic radii set to "csd" values;
  * "alv": Alvarez version;
  * "chap": Chernyshov&Ananyev&Pidko (this work) version.
* *--lab1*: name of CSV file column containing labels of the contact's first atom, "LAB1" by default;
* *--lab2*: name of CSV file column containing labels of the contact's second atom, "LAB2" by default;
* *--dist*: name of CSV file column containing contact distances, "DIST1" by default;

**Please make sure, that *--lab1*, *--lab2* and *--dist* values correspond to those in CSV file.**

* *--norm*: type of X–H bonds normalization, "csd" by default. **Please make sure, that the same normalization scheme was used for a ConQuest search.** Available values are:
  * "csd": C–H, N–H and O–H bonds are normalized to 1.089 Å, 1.015 Å and 0.993 Å, correspondingly;
  * "no": X–H bonds are not normalized;
  * path to the file, each line of those contains space separated element symbol and the length of corresponding X–H bond in angstroms, e.g. "C 1.09".
* *--tol*: minimal possible distance between different atoms, default 0.005 Å.
* *-V* or *--volume*: maximal allowed cell volume, Å<sup>3</sup>. Available values are:
  * positive numeric value: in this case the script does not calculate contact shielding for crystals with crystallographic cell volume more than specified. It can be useful if there are a lot of crystals with large crystallographic cells (V > 10 000 Å<sup>3</sup>) which are treated slowly.
  * unspecified: no filtering by volume is applied.

