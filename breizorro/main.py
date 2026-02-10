import os

import click
from omegaconf import OmegaConf
from scabha.schema_utils import clickify_parameters

from breizorro.breizorro import main

schemas = OmegaConf.load(os.path.join(os.path.dirname(__file__), "breizorro.yaml"))


@click.command("breizorro")
@clickify_parameters(schemas.cabs.get("breizorro"))
def driver(**kw):
    main(**kw)
