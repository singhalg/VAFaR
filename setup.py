#-------------------------------------------------------------------------------
# Name:        setup.py
# Purpose:
#
# Author:      Gaurav Singhal
#
# Created:     12/06/2013
# Copyright:   (c) Gaurav Singhal 2013
# Licence:     This work is licensed under the Creative Commons Attribution-NonCommercial 3.0 Unported License.
#              To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/3.0/.
#-------------------------------------------------------------------------------

from setuptools import setup

setup(name='VAFaR',
      version='1',
      description='Variant Analysis, Filtration and Representation',
      url='http://github.com/singhalg/VAFaR',
      author='Gaurav Singhal',
      author_email='gaurav.singhal.01@gmail.com',
      license='Creative Commons Attribution-NonCommercial 3.0 Unported License',
      packages='vafar',
      install_requires=['set','pickle', 'copy', 'numpy', 'subprocess', 'optparse', 'sys', 'random', 'math'])

def main():
    pass

if __name__ == '__main__':
    main()
