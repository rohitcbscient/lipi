import numpy as np
import matplotlib.pyplot as plt
def get_tt(r: int):
    """_summary_

    Args:
        r: _description_
    """
    print(r*100)
    return 'r'
a=get_tt(3)

"""
[project]
name = "deepcanary"
description = Real-time screen watermarking to determine source of leaks.
authors = [
    { name = "Lukas Gehrig", email = "lukas.gehrig@fhnw.ch" },
    { name = "Cédric Huwyler", email = "cedric.huwyler@fhnw.ch" },
]
readme = "README.md"
requires-python = ">=3.12"
dynamic = ["version"]

dependencies = [
    "bitstring >= 4.2.3",
    "bchlib @ git+https://github.com/m472/python-bchlib@0.14.1",  # install updated bchlib from PyPi once `bch.decode` has been fixed in 1.0.0
    "matplotlib >= 3.9.2",
    "numpy >= 1.26.4",
    "opencv-python >= 4.10.0",
    "pandas >= 2.2.2",
    "pyserde[all] >= 0.20.1",
    "requests >= 2.32.3",
    "scipy >= 1.14.0",
    "torch >= 2.6.0",  # CUDA 12.4 without index-url [10-2024]
    "torchmetrics >= 1.4.1",
    "torchvision >= 0.19.0",
]

[project.optional-dependencies]

    train = [
        "kornia >= 0.7.3",
        "lpips >= 0.1.4",
        "wandb >= 0.17.7",
    ]

	test = [
        "pytest == 8.1.1",
        "pytest-xdist >= 3.6.1",
    ]
    
    dev = [
		"argcomplete == 3.5.0",
        "mypy == 1.15.0",
        "pre-commit >= 3.8.0",
		"rich == 13.7.1",
        "ruff == 0.6.1",
        "types-requests >= 2.32.0",
        "types-setuptools >= 71.1.0",
    ]

    jupyter = [
        "ipykernel >= 6.29.5",
        "ipylab >= 1.0.0",
        "ipywidgets >= 8.1.3",
        "jupyterlab >= 4.2.4",
    ]

[build-system]
requires = ["setuptools>=56.0", "wheel", "versioneer[toml]"]
build-backend = "setuptools.build_meta"

[tool.setuptools.dynamic]
version = { attr = "deepcanary.__version__" }

[tool.mypy]
ignore_missing_imports = true
warn_unused_configs = true
warn_redundant_casts = true
warn_unused_ignores = true
no_implicit_optional = true
strict_equality = true
extra_checks = true
check_untyped_defs = true
disallow_subclassing_any = true
disallow_untyped_decorators = true
disallow_any_generics = true
disallow_untyped_calls = true
disallow_incomplete_defs = true
disallow_untyped_defs = true
exclude = [
    '^.*test_.*\.py$',
]

[tool.pyright]
reportImplicitStringConcatenation = true
reportMissingParameterType = true
reportMissingSuperCall = true
reportUnnecessaryTypeIgnoreComment = true
reportUnknownParameterType = true

[tool.versioneer]
VCS = "git"
style = "pep440"
versionfile_source = "src/deepcanary/_version.py"
versionfile_build = "deepcanary/_version.py"
tag_prefix = "v"
parentdir_prefix = "deepcanary-"

[tool.pytest.ini_options]
testpaths = ["test"]
markers = ["long", "training"]

[tool.ruff]
required-version = ">=0.6.1"
line-length = 88
exclude = ["versioneer.py", "src/deepcanary/_version.py"]
extend-include = ["*.ipynb"]
src = ["src", "test"]

[tool.ruff.lint]
extend-select = ["E501"]
logger-objects = ["deepcanary.utils.logger"]
exclude = ["__init__.py", "test/*"]
select = [  # see https://docs.astral.sh/ruff/rules/
    "F",  # Pyflakes
    "E",  # pycodestyle/Error
    "W",  # pycodestyle/Warning
    "I",  # isort
    "N",  # pep8-naming
    "D",  # pydocstyle
    "UP",  # pyupgrade
    "ANN",  # flake8-annotations
    "B",  # flake8-bugbear
    "C4",  # flake8-comprehensions
    "EM",  # flake8-errmsg
    "FA",  # flake8-future-annotations
    "ICN",  # flake8-import-conventions
    "LOG",  # flake8-logging
    "INP",  # flake8-no-pep420
    "PIE",  # flake8-pie
    "PYI",  # flake8-pyi
    "PT",  # flake8-pytest-style
    "Q",  # flake8-quotes
    "RET501",  # flake8-return unnecessary-return-none
    "RET502",  # flake8-return implicit-return-value
    "RET504",  # flake8-return unnecessary-assign
    "SLOT",  # flake8-slots
    "SIM",  # flake8-simplify
    "TID",  # flake8-tidy-imports
    "TCH",  # flake8-type-checking
    "ARG",  # flake8-unused-arguments
    "ERA",  # eradicate
    "PD",  # pandas-vet
    "PGH",  # pygrep-hooks
    "PLC",  # Pylint/Convention
    "PLE",  # Pylint/Error
    "PLW",  # Pylint/Warning
    "TRY",  # tryceratops
    "NPY",  # NumPy-specific rules
    "PERF",  # Perflint
    "RUF",  # Ruff-specific rules
]
ignore = [
    "D105",  # docstring in magic-method
    "UP040",  # PEP 695 not supported yet by mypy (04-2024)
    "ANN101",  # missing-type-self is deprecated
    "ANN102",  # missing-type-cls is deprecated
]

[tool.ruff.format]
docstring-code-format = true

[tool.ruff.lint.flake8-annotations]
allow-star-arg-any = true

[tool.ruff.lint.flake8-tidy-imports]
ban-relative-imports = "all"

[tool.ruff.lint.pycodestyle]
ignore-overlong-task-comments = true
max-line-length = 99

[tool.ruff.lint.pydocstyle]
convention = "google"
ignore-decorators = ["typing.overload"]
"""