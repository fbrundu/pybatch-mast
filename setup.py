#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'boto3>=1.14.20',
    'pandas>=1.1.2',
    'scanpy>=1.6.0',
    'XlsxWriter==1.2.9',
]

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Francesco G. Brundu",
    author_email='francesco.brundu@gmail.com',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Differential expression analysis with MAST on AWS Batch (python3 interface).",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pybatch_mast',
    name='pybatch_mast',
    packages=find_packages(include=['pybatch_mast', 'pybatch_mast.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/fbrundu/pybatch_mast',
    version='0.1.0',
    zip_safe=False,
)
