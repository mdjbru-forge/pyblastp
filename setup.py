# Written usingn resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages

setup(name = "pyblastp",
      version = "0.0.1",
      py_modules = ["pyblastp", "pyblastpScripts"],
      install_requires = ["biopython"],
      entry_points =  {
          "console_scripts" : [
              "pyblastp=pyblastpScripts:main"
          ]
      }
)
