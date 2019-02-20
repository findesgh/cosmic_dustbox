from setuptools import setup

with open('README.md') as fh:
    long_desc = fh.read()
name = 'cosmic_dustbox'

setup(
    name=name,
    version='0.1.dev0',
    author_email='findessp@yandex.ru',
    description='Toolbox for cosmic dust.',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    url='https://github.com/findesgh/'+name,
    license='GPL-3.0',
    packages=[name],
    install_requires=[
        'numpy',
        'astropy',
        'scipy'
    ],
    tests_require=[
        'pytest'
    ],
    zip_safe=False,
)
