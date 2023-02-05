<h1 align="center">py2sambvca</h1> 
<h3 align="center">Simple thin client to interface python scripts with SambVca catalytic pocket Fortran calculator.</h3>

<p align="center">
  <img alt="py2sambvca logo" src="https://github.com/JacksonBurns/py2sambvca/blob/main/py2sambvca_logo.png">
</p>

<p align="center">
  <img alt="GitHub Repo Stars" src="https://img.shields.io/github/stars/JacksonBurns/py2sambvca?style=social">
  <img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/py2sambvca">
  <img alt="Total Downloads" src="https://static.pepy.tech/personalized-badge/py2sambvca?period=total&units=international_system&left_color=grey&right_color=blue&left_text=Downloads">
  <img alt="PyPI" src="https://img.shields.io/pypi/v/py2sambvca">
  <img alt="commits since" src="https://img.shields.io/github/commits-since/JacksonBurns/py2sambvca/latest.svg">
  <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/py2sambvca">
  <img alt="PyPI - License" src="https://img.shields.io/pypi/l/py2sambvca">
</p>

## Installation
`py2sambvca` is available on PyPI and can be installed like so:
```python
pip install py2sambvca
```

`py2sambvca` has __zero__ external depdencies.

### Downloading and Compiling `Sambvca`
`py2sambvca` can read and write input and output files for `Sambvca` without the actual program in place, but in order to run input files you must have an executable `sambvca21.exe` (or similar) somewhere on your machine.

You can download the source code [on the `Sambvca` webserver](https://www.molnac.unisa.it/OMtools/sambvca2.1/download/download.html) and compile it using [`gfortran`](https://gcc.gnu.org/wiki/GFortranBinaries).

By default, `py2sambvca` expects the executable to be present in the `cwd` and named `sambvca21.exe` on Windows or `sambvca21.x` on Unix-based systems. optionally, the filepath to your executable can be specified as shown below.

## Usage
After installation, `py2sambvca` can be added to a Python script via `import` and instantiated:
```python
from py2sambvca import p2s

nhc_p2s = p2s(
    "test/data/nhc.xyz",
    [22],
    [5],
    [1],
    path_to_sambvcax="sambvca21.exe",
)
```

The required input parameters are shown below:
 - `xyz_filepath` (str): Location of .xyz molecular coordinates file for writing input data
 - `sphere_center_atom_ids` (list): ID of atoms defining the sphere center
 - `z_ax_atom_ids` (list): ID of atoms for z-axis
 - `xz_plane_atoms_ids` (list): ID of atoms for xz-plane

The following parameters are optional and will be filled with default values if not specified:
 - `atoms_to_delete_ids` (list): ID of atoms to be deleted (default None)
 - `sphere_radius` (float): Sphere radius in Angstrom (default 3.5)
 - `displacement` (float): Displacement of oriented molecule from sphere center in Angstrom (default 0.0)
 - `mesh_size` (float): Mesh size for numerical integration (default 0.10)
 - `remove_H` (int): 0/1 Do not remove/remove H atoms from Vbur calculation (default 1)
 - `orient_z` (int): 0/1 Molecule oriented along negative/positive Z-axis (default 1)
 - `write_surf_files` (int): 0/1 Do not write/write files for top and bottom surfaces (default 1)
 - `path_to_sambvcax` (str): Path to the SambVca executable. Only needed to use py2sambvca.calc()(default "sambvca.exe")
 - `working_dir` (path): Path to the working directory where the output and input files are generated (default os.getcwd())
 - `verbose` (int): 0 for no output, 1 for some output, 2 for the most output (default 1)


From here, running can be done stepwise or with a single function:
```python
nhc_p2s.run()
# equivalent to
nhc_p2s.write_input()
nhc_p2s.calc()
nhc_p2s.parse_output()
nhc_p2s.clean_files()
```

All values for the total complex, quadrants, and octants are available through getters:

| Total Values | Quadrant Values | Octant Values |
| --- | --- | --- |
| `get_free_volume()` | `get_quadrant_free_volume()` | `get_octant_free_volume()` |
| `get_buried_volume()` | `get_quadrant_buried_volume()` | `get_octant_buried_volume()` |
| `get_exact_volume()` | _not available_ | _not available_ |
| `get_total_volume()` | `get_quadrant_total_volume()` | `get_octant_total_volume()` |
| `get_percent_buried_volume()` | `get_quadrant_percent_buried_volume()` | `get_octant_percent_buried_volume()` |
| `get_percent_free_volume()` | `get_quadrant_percent_free_volume()` | `get_octant_percent_free_volume()` |
| `get_percent_total_volume()` | _not available_ | _not available_ |

Results can also be accessed through a general getter method: `get()`, `get_quadrant_result()`, and `get_octant_result()`.

All results can also be directly accessed through dictionaries, returned from a call to `run()` or `parse_output()` and available through `p2s.total_results`, `p2s.quadrant_results`, and `p2s.octant_results`.

In case there is something else you are looking for, you can use a general purpose `get_regex()` function to return the line containing a pattern.

### Examples
Here are a couple repositories using `py2sambvca` as a Python package or extending its source code, check them out:
 - ~~[Metal-organic framework stability analysis by Hiu Ki](https://github.com/hiukiwong/mof-stability-ml)~~
 - [MOF Stability ML by Ruihan Wang](https://github.com/ruihwang/mof-stability-ml)

### See Also
 - Kjell Jorner's [morfeus](https://github.com/kjelljorner/morfeus) package re-implements the original buried volume algorithm directly in Python

## License
`py2sambvca` is available under the GNU GPLv3 in accordance with the base Fortran code which is available under the same license and can be retreieved here: https://www.molnac.unisa.it/OMtools/sambvca2.1/download/download.html

The original fortran program (`sambvca21.f`) is also included in the `test` directory for testing purposes. It is still under the same terms of the GNU license:
 - This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
 - This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 - The results obtained from using the source code shall be used for scientific purposes only, excluding industrial or commercial purposes. To use the SambVca suite for industrial or commercial purposes, contact lcavallo|@|unisa.it.
 - Proper acknowledgement shall be made to the author of the source code in publications resulting from the use of it in its original form or modified.
 - The results from using the source code are provided "AS IS" without warranty of any kind.

## Citation
Please cite the `SambVca` underlying Fortran tool according to the guidelines on the buried volume webserver: [https://www.molnac.unisa.it/OMtools/sambvca2.1/help/help.html](https://www.molnac.unisa.it/OMtools/sambvca2.1/help/help.html)

`py2sambvca` has been uploaded to Figshare and may be cited as: Burns, J. figshare. 2020, DOI:[10.6084/m9.figshare.12846707](https://figshare.com/articles/software/py2sambvca/12846707)
