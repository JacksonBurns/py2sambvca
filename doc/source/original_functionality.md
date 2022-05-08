# Usage

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

# Examples

Here are a couple repositories using `py2sambvca` as a Python package or extending its source code, check them out:

- [Metal-organic framework stability analysis by Hiu Ki](https://github.com/hiukiwong/mof-stability-ml)
- [Phosphine ligand parameterization for Machine Learning by Kjell Jorner](https://github.com/kjelljorner/morfeus)
