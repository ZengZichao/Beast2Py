"""Beast2Py: A Python framework for reproducible divergence time
estimation with automated calibration prior specification, validation, and diagnostics.

Beast2Py provides two calling modes:

1. **CLI mode**: Use the ``beast2py`` command-line tool::

       beast2py generate --config config.yaml --output output.xml
       beast2py diagnose --config config.yaml --report report.html

2. **API mode**: Import Beast2Py as a Python library::

       from beast2py.api import Beast2Py

       b2p = Beast2Py()
       xml = b2p.generate_xml("config.yaml", output="output.xml")
       report = b2p.diagnose("config.yaml", report_path="report.html")

Both modes provide identical functionality. The API is recommended for
programmatic use, pipeline integration, and Jupyter notebook workflows.
"""

__version__ = "0.1.1"
__author__ = "曾子超 (Zichao Zeng)"
__email__ = "zengzichao@sjtu.edu.cn"
__orcid__ = "0000-0001-6553-970X"
__license__ = "MIT"

# Convenience imports for common API usage
from .api import (
    Beast2Py,
    generate_xml,
    quick_generate,
    diagnose,
    validate_xml,
    generate_fingerprint,
    generate_methods,
    generate_pipeline,
    list_models,
    parse_config,
    read_sequence,
)

# Re-export key types for type hints
from .models import (
    BEASTConfig,
    CalibrationMethod,
    CalibrationPoint,
    ClockModelConfig,
    ClockModelType,
    DataType,
    DistributionConfig,
    InitTreeType,
    InitializationConfig,
    LoggerConfig,
    MCMCConfig,
    MCMCType,
    Partition,
    Provenance,
    RealParameter,
    SiteModelConfig,
    TipDatesConfig,
    TreePriorType,
)

__all__ = [
    # API class
    "Beast2Py",
    # Module-level API functions
    "generate_xml",
    "quick_generate",
    "diagnose",
    "validate_xml",
    "generate_fingerprint",
    "generate_methods",
    "generate_pipeline",
    "list_models",
    "parse_config",
    "read_sequence",
    # Configuration types
    "BEASTConfig",
    "CalibrationMethod",
    "CalibrationPoint",
    "ClockModelConfig",
    "ClockModelType",
    "DataType",
    "DistributionConfig",
    "InitTreeType",
    "InitializationConfig",
    "LoggerConfig",
    "MCMCConfig",
    "MCMCType",
    "Partition",
    "Provenance",
    "RealParameter",
    "SiteModelConfig",
    "TipDatesConfig",
    "TreePriorType",
]
