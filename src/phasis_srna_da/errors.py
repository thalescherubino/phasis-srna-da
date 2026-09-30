"""Domain-specific exceptions with messages intended for command-line users."""


class ValidationError(ValueError):
    """Raised when an input violates the documented analysis contract."""


class ExternalToolError(RuntimeError):
    """Raised when Bowtie or Bowtie-build fails."""
