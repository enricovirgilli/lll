class plot2d(object):
    """Class defining material properties useful for crystallography."""
    def __init__(self, *args, **keyargs):
        self.files=args
        for key in keyargs.keys(): self.__setattr__(key, keyargs[key])

    def __str__(self):
        """Class representation"""
        return "%s, fwhm=%s arcmin, microblock thickens=%s" % (Booklet.name[self.Z], self.fwhm, self.microthick)

