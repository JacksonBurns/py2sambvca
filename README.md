# py2sambvca
![GitHub Repo stars](https://img.shields.io/github/stars/JacksonBurns/py2sambvca?style=social)
![PyPI - Downloads](https://img.shields.io/pypi/dm/py2sambvca)
![PyPI](https://img.shields.io/pypi/v/py2sambvca)
![PyPI - License](https://img.shields.io/pypi/l/py2sambvca)

 Simple thin client to interface python scripts with SambVca catalytic pocket Fortran calculator.

## Installation
`py2sambvca` is available on PyPi and can be installed like so:
```python
pip install py2sambvca
```

`py2sambvca` has __zero__ external depdencies.

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

Total Values:
 - `get_free_volume()`
 - `get_buried_volume()`
 - `get_exact_volume()`
 - `get_total_volume()`
 - `get_percent_buried_volume()`
 - `get_percent_free_volume()`
 - `get_percent_total_volume()`

Quadrant Values:
 - `get_quadrant_free_volume()`
 - `get_quadrant_buried_volume()`
 - `get_quadrant_total_volume()`
 - `get_quadrant_percent_buried_volume()`
 - `get_quadrant_percent_free_volume()`

Octant Values:
 - `get_octant_free_volume()`
 - `get_octant_buried_volume()`
 - `get_octant_total_volume()`
 - `get_octant_percent_buried_volume()`
 - `get_octant_percent_free_volume()`

Results can also be accessed through a general getter method: `get()`, `get_quadrant_result()`, and `get_octant_result()`.

All results can also be directly accessed through dictionaries, returned from a call to `run()` or `parse_output()` and availabel through `p2s.total_results`, `p2s.quadrant_results`, and `p2s.octant_results`.

In case there is something else you are looking for, you can use a general purpose `get_regex()` function to return the line containing a pattern.

### Examples
Here are a couple repositories using `py2sambvca` as a Python package or extending its source code, check them out:
 - ~~[Metal-organic framework stability analysis by Hiu Ki](https://github.com/hiukiwong/mof-stability-ml)~~
 - [Phosphine ligand parameterization for Machine Learning by Kjell Jorner](https://github.com/kjelljorner/morfeus)
 - [MOF Stability ML by Ruihan Wang](https://github.com/ruihwang/mof-stability-ml)

## License
`py2sambvca` is available under the GNU GPLv3 in accordance with the base Fortran code which is available under the same license and can be retreieved here: https://www.molnac.unisa.it/OMtools/sambvca2.1/download/download.html

The original fortran program (`sambvca21.f`) is also included in the `test` directory for testing purposes. It is still under the same terms of the GNU license:
 - This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.
 - This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 - The results obtained from using the source code shall be used for scientific purposes only, excluding industrial or commercial purposes. To use the SambVca suite for industrial or commercial purposes, contact lcavallo|@|unisa.it.
 - Proper acknowledgement shall be made to the author of the source code in publications resulting from the use of it in its original form or modified.
 - The results from using the source code are provided "AS IS" without warranty of any kind.

## Citation
Please cite the `SambVca` base fortran tool as: Falivene, L. et al. Nat. Chem. 2019, DOI:10.1038/s41557-019-0319-5 

`py2sambvca` has been uploaded to Figshare and may be cited as: Burns, J. figshare. 2020, DOI:10.6084/m9.figshare.12846707
