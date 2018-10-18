import os

def rv_format(fn, kind='i', **kwds):
    """Read the model input template file and fill the given arguments."""

    if None in kwds.values():
        raise ValueError("Some parameters are not properly set.")

    if kind not in ['chipt']:
        raise ValueError('`kind` must be one of c,h,i,p,t: {}'.format(kind))

    with open(fn + '.rv' + kind) as f:
        txt = f.read()

    return txt.format(**kwds)


def create_subdir(name, outpath, **kwds):
    """Create a subdirectory and write the RAVEN configuration files. """
    from raven.models import raven_templates
    inpath = raven_templates[name]

    # Create subdirectory
    os.mkdir(os.path.join(outpath, 'model'))
    os.mkdir(os.path.join(outpath, 'output'))

    for kind in 'chipt':
        txt = rv_format(name, kind, **kwds)
        with open(os.path.join(outpath, model, name + '.rv' + kind), 'w') as f:
            f.write(txt)

    os.symlink(os.path.join(..., 'raven_rev.exe'))

