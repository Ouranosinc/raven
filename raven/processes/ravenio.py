def rv_format(fn, kind='i', **kwds):
    """Read the model input template file and fill the given arguments."""

    if None in kwds.values():
        raise ValueError("Some parameters are not properly set.")

    if kind not in ['chipt']:
        raise ValueError('`kind` must be one of c,h,i,p,t: {}'.format(kind))

    with open(fn + '.rv' + kind) as f:
        txt = f.read()

    return txt.format(**kwds)
