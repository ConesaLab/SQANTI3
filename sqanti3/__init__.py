from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass

__author__ = __email__ = ("fraparp1@upv.edu.es", "pedsalga@upv.edu.es")
