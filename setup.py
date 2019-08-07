from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='geneLift',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Gene model transfer from closely related reference genomes using cDNA alignments",
    author="Srividya",
    author_email='srividya.ramki@gmail.com',
    url='https://github.com/srividya22/geneLift',
    packages=['geneLift'],
    entry_points={
        'console_scripts': [
            'geneLift=geneLift.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='geneLift',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ]
)
