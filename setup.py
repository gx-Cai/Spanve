from setuptools import setup, find_packages
import os

# Get the absolute path to the README
here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'ReadMe.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='Spanve',
    version='0.1.6',
    description='Spatial Neighborhood Variably Expressed Genes (Spanve) is a method for detecting spatial variably expressed genes in spatial transcriptomics data.',
    url='https://github.com/gx-Cai/Spanve',
    author='gx.Cai',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy',
        'scanpy',
        'scipy>=1.8',
        'joblib',
        'scikit-learn',
        'statsmodels',
        'tqdm',
        'click'
    ],
    entry_points={
        'console_scripts': [
            'Spanve = Spanve.Spanve_cli:main'
        ]
    },
    long_description=long_description,
    long_description_content_type='text/markdown',
)

