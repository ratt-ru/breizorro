import os
from setuptools import setup, find_packages

pkg = 'breizorro'
__version__ = "0.0.1"
build_root = os.path.dirname(__file__)

def readme():
    """Get readme content for package long description"""
    with open(os.path.join(build_root, 'README.rst')) as f:
        return f.read()

def requirements():
    """Get package requirements"""
    with open(os.path.join(build_root, 'requirements.txt')) as f:
        return [pname.strip() for pname in f.readlines()]

setup(name=pkg,
      version=__version__,
      description="Creates a binary mask given a FITS image",
      long_description=readme(),
      author="Ian Heywood & RATT",
      author_email="ian.heywood@physics.ox.ac.uk",
      packages=find_packages(),
      url="https://github.com/ratt-ru/breizorro",
      license="GNU GPL 3",
      classifiers=["Development Status :: 4 - Beta",
                   "Intended Audience :: Developers",
                   "Programming Language :: Python :: 3.6",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   "Topic :: Software Development :: Libraries :: Python Modules"],
      keywords="fits dataset models mask",
      platforms=["OS Independent"],
      install_requires=requirements(),
      python_requires='>=3.6',
      include_package_data=True,
      scripts=['breizorro/bin/breizorro'])
