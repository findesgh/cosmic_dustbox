from setuptools import setup

with open('README.md') as fh:
    long_desc = fh.read()

setup(
    name='cosmic_dustbox',
    version='0.1.dev',
    author_email='findessp@yandex.ru',
    description='Toolbox for cosmic dust.',
    long_description=long_desc,
    long_description_content_type='text/markdown',
    url='https://github.com/findesgh/cosmic_dustbox',
    license='GPL-3.0',
    packages=['cosmic_dustbox'],
    install_requires=[
        'numpy'
    ],
    zip_safe=False,
)
