from py2sambvca.radii_tables.default import default_radii_table
from py2sambvca.radii_tables.vanDerWaals import vdw_radii_table
from py2sambvca.radii_tables.format_table import format_radii_table

table_lookup = {
    "default": default_radii_table,
    "vdw": vdw_radii_table,
}

__all__ = ["format_radii_table", "table_lookup"]
