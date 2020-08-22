# py2sambvca
 Simple thin client to interface python scripts with SambVca catalytic pocket Fortran calculator.

## Installation
`py2sambvca` is available on PyPi and can be installed like so:
```python
pip install py2sambvca
```

## Usage
After installation, `py2sambvca` can be added to a Python script via `import` and instantiated:
```python
# import the class and give it a simpler alias
from py2sambvca.py2sambvca import py2sambvca as p2s

buried_vol = p2s(r'myxyzfiles\ligand_4.xyz',...)
buried_vol.write_input()
buried_vol.calc()

# retrieve the buried volume and assign it
ligand_4_buried_volume = buried_vol.get_buried_vol()
...
# clean up input and output files
buried_vol.clean_files()
```

## License
`py2sambvca` is available under the GNU GPLv3 in accordance with the base Fortran code which is available under the same license and can be retreieved here: https://www.molnac.unisa.it/OMtools/sambvca2.1/download/download.html
