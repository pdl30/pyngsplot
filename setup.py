import os
from setuptools import setup, find_packages

setup(name='pyngsplot',
      version='0.0.1',
      description='pyngsplot is a Python module to analyze plot heatmaps and average profiles of NGS aligned data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=find_packages(),
      package_data={"pyngsplot":['data/*']},
      scripts=['scripts/pyngsplot.py'],
      install_requires=['pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pyngsplot is a Python module to analyze plot heatmaps and average profiles of NGS aligned data

 Contact
=============

If you have any questions or comments about pyngsplot, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
