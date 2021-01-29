from setuptools import setup

print('This 3.6+ script requires rdkit, numpy and optionally pymol and jupyter, which are best installed with conda.')

setup(
    name='Fragmenstein',
    version='0.5',
    packages=['fragmenstein'],
    install_requires=['numpy', 'rdkit'],
    extras_require={'minimization': ['pyrosetta','rdkit_to_params'],
                      'pse-export': ['pymol'],
                      'jupyter': ['jupyter']},
    url='https://github.com/matteoferla/Fragmenstein',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Scaffold hopping between bound compounds by stitching them together like a reanimated corpse',
    entry_points={
        'console_scripts': ['fragmenstein=cli/fragmenstein.py:main'],
    }
)
