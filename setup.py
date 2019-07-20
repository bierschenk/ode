from setuptools import setup

with open('README.rst', 'r') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst', 'r') as history_file:
    history = history_file.read()

setup(
    name='ode',
    description="Ordinary differential equation solver (numeric integration)",
    long_description=readme + '\n\n' + history,
    keywords='ode',
    version='0.3.0',
    license="MIT license",
    packages=['ode'],
    include_package_data=True,
    test_suite='tests',
    author="bierschenk",
    author_email='bierschenk.devel@gmail.com',
    url='https://github.com/bierschenk/ode',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering'
    ],
    install_requires=[],
)
