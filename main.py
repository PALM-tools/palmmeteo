from argparse import ArgumentParser
from core import run


argp = ArgumentParser()
argp.add_argument('-c', '--config', default='config.yaml')

argv = argp.parse_args()

run(argv)
