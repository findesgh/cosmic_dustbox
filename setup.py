from setuptools import setup
import metadata as md


setup(
    name=md.name,
    version=md.version,
    author=md.authors,
    author_email=md.author_email,
    description=md.description,
    long_description=md.long_desc,
    long_description_content_type='text/markdown',
    url=md.url,
    license=md.license,
    packages=[md.name],
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
