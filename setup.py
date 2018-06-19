from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys

requirements = [line.strip() for line in open('requirements.txt').readlines()]

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest


setup(name='screening_template',
      version='0',
      description='',
      url='https://github.com/ahy3nz/screening_template',
      author='Alexander Yang',
      author_email='alexander.h.yang@vanderbilt.edu',
      license='MIT',
      packages=['screening_template'],
      install_requires=requirements,
      zip_safe=False,
      test_suite='tests',
      cmdclass={'test': PyTest},
      #extras_require={'utils': ['pytest']},
)
