from collections import OrderedDict


# TODO: Implement section parser
def parse_configuration(fn):
    """Parse Raven configuration file.

    Returns a dictionary keyed by parameter name."""
    import re

    main_param = re.compile(r"^:(\w+)\s+([^#]*)")
    # sub_param = re.compile(r"^  :(\w+)\s+([^#]*)")
    out = OrderedDict()
    # cat = None
    with open(str(fn)) as f:
        for line in f.readlines():
            match = main_param.search(line)
            if not match:
                continue

            key, value = match.groups()
            if value:
                values = value.split()
                out[key] = values[0] if len(values) == 1 else values
            else:
                if "List" in key:
                    pass
                elif "Classes" in key:
                    pass
                elif "Profiles" in key:
                    pass
                else:
                    out[key] = True

    return out
