import sys

__all__ = ['die', 'warn', 'log', 'verbose', 'configure_log']


def die(s, *args, **kwargs):
    """Write message to error output and exit with status 1."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')
    sys.exit(1)


def warn(s, *args, **kwargs):
    """Write message to error output."""

    if args or kwargs:
        error_output(s.format(*args, **kwargs) + '\n')
    else:
        error_output(s + '\n')


def log_on(s, *args, **kwargs):
    """Write logging or debugging message to standard output (logging is enabled)."""

    if args or kwargs:
        log_output(s.format(*args, **kwargs) + '\n')
    else:
        log_output(s + '\n')


def log_off(s, *args, **kwargs):
    """Do nothing (logging is disabled)."""

    pass


log_output = sys.stdout.write
error_output = sys.stderr.write
log = log_on
verbose = log_off

def configure_log(cfg):
    global log, verbose

    log = log_on if cfg.verbosity >= 1 else log_off
    verbose = log_on if cfg.verbosity >= 2 else log_off
