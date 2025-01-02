from setuptools import setup, find_packages

setup(
    name='monet',
    version='0.1',
    author='Nimrod Rappoport',
    url='https://github.com/Shamir-Lab/MONET',
    description='Python implementation of the MONET algorithm for multi-omic module detection',
    license='GPL-3.0',
    python_requires='>=3.0.0',
    install_requires=[
        'numpy',
        'pandas',
        'networkx',
        'scipy'
    ],
    packages=['monet'],  # Explicitly specify the package to include
    # include_package_data=True,  # Include non-Python files if needed
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
)

