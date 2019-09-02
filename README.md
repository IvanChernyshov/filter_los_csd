# filter_los_csd

This is the manual for the script which filters line-of-sight contacts from the results of ConQuest search over the Cambridge Structure Database (CSD) entries. The script is available via [GitHub](https://github.com/IvanChernyshov/filter_los_csd). For the underlying scientific context see the corresponding work.

If you use the script in your work, please cite generously!

This manual consists of two parts:

* *Installation* containing script requirements and example of conda installation;
* *Syntax* of the script.

The example of input and output files for the script is provided in *example* directory.

## Installation

*filter_los_csd* requires installation of:

1. Python 3;
2. [numpy](https://numpy.org/);
3. [pandas](https://pandas.pydata.org/);
4. [PyCifRW](https://www.iucr.org/resources/cif/software/pycifrw). 

All this libraries can be installed on Windows, Linux and macOS. Here we provide installation using cross-platform software package manager Conda.

The easiest way to get Conda is having it installed as part of the [Anaconda Python distribution](https://www.anaconda.com/distribution/). A possible (but a bit more complex to use) alternative is provided with the smaller and more self-contained [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

After conda installation, create new environment and install necessary packages:

```
> conda create -n los python=3 numpy pandas
> conda install -n los -c conda-forge pycifrw
```

Finally, the new environment must be activated so that the corresponding python interpreter becomes available in the same shell:

```
> source activate los
```

or, in Windows:

```
> activate los
```

For more details see the [Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).

## Syntax

### Input and Output

The script takes two necessary parameters as input:

* *path_csv*: path to CSV file, containing info on parameters of contacts in question.
  **Please note, that CSV file must contain labels of atoms forming the contact and the contact distance!**
* *path_cif*: path to multiple CIF file, containing CSD entries from CSV file.
  **Please note, that CIF file must contain info on bonds.** Make sure that the corresponding checkbox was selected before downloading CIF.

Thus, the easiest command is:

```
> python filter_csd_los.py test.csv test.cif
```

As output the script creates **{csv_name}_los.csv** file containing the same info as the original file with three additional columns:

* "LOS": has three possible values:
  * "+": the corresponding contact is line-of-sight;
  * "–": the corresponding contact is not line-of-sight;
  * "?": the error occurred during the calculation. The only type of errors caught during the testing is inability to find a contact with given atomic labels and distance. The main source of these errors are disorder issues, thus tuning *--tol* parameter can solve the problem.
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

