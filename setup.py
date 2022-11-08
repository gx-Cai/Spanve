from setuptools import setup, find_packages

setup(
    name='Spanve',
    version='0.1.0',
    description='Spatial Neighbourhood Variably Expressed (Spanve) is a method for detecting spatially expressed genes in spatial transcriptomics data.',
    url='github.com/gx-Cai/Spanve',
    author='gx.Cai',
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy<=1.21',
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
    }
)

