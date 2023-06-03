# this radius table just returns to the plain old van Der Waals radii,
# which are the radii reported in the original Cavallo paper divided by
# 1.17
from py2sambvca.radii_tables.default import default_radii_table

vdw_radii_table = {
    element: radius / 1.17 for element, radius in default_radii_table.items()
}
