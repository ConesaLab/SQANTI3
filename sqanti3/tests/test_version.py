import unittest
import logging

try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

logging.basicConfig(level=logging.CRITICAL)


class TestVersion(unittest.TestCase):
    def test_version(self):
        print(metadata.version('SQANTI3'))


if __name__ == "__main__":
    unittest.main()
