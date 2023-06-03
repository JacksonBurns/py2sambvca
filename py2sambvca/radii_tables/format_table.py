def format_radii_table(radii_table: dict):
    """Formats a dictionary of atomic symbol: radii to a list of strings for sambvca

    Args:
        radii_table (dict): Mapping of atomic symbols to radii, i.e. {'H':1.0,'LI':1.1}
    """
    output_list = []
    for element, radius in radii_table.items():
        # sambvca expects each row in the radii table to look like this:
        # "C       1.11" with _exactly_ six spaces for elements with two letter
        # abbreviations and _exactly_ seven for elements with one letter
        # abbreviations, and then three digits for the radius
        output_list.append(
            element
            + (" " if len(element) == 1 else "")
            + "      "
            + str(round(radius, 2))
            + "\n"
        )
    return output_list
