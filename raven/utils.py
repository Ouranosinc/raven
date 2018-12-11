from re import search


def address_append(address):
    zipped = search(r'(\.zip)', address)
    tarred = search(r'(\.tar)', address)

    try:
        if zipped:
            return 'zip://{}'.format(address)
        elif tarred:
            return 'tar://{}'.format(address)
        else:
            return address
    except Exception as e:
        msg = 'Failed to prefix or parse URL: {}'.format(e)
        raise Exception(msg)
