#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import os
from pathlib import Path
from typing import List

import pytest
import nbmake


class KernelSpecMismatch(ValueError):
    def __init__(self, requested: str, specified: str):
        self.requested = requested
        self.specified = specified
        super().__init__(
            f"Kernelspec with name {self.requested!r} requested, but notebook specifies {self.specified!r}"
        )


class KernelSpec(pytest.Item):
    def __init__(self, requested: str = "", **kwargs):
        super().__init__(**kwargs)
        self.requested = requested

    def setup(self):
        import nbformat

        self._nb = nbformat.read(self.path, as_version=4)
        self._metadata = self._nb.metadata
        self._kernelspec = self._metadata.kernelspec

    def runtest(self):
        if self._kernelspec.name != self.requested:
            raise KernelSpecMismatch(
                requested=self.requested,
                specified=self._kernelspec,
            )

    def repr_failure(self, excinfo):
        exc = excinfo.value
        if isinstance(exc, KernelSpecMismatch):
            return str(exc)
        return super().repr_failure(excinfo)

    def reportinfo(self):
        return (
            str(self.path),
            None,
            self.name,
        )


class NotebookMetadata(pytest.File):
    def collect(self):
        requested_kernel = self.parent.config.option.nbmake_kernel
        if requested_kernel:
            yield KernelSpec.from_parent(
                self,
                name=f"{self.name}:kernelspec",
                path=self.path,
                requested=requested_kernel,
            )


def pytest_collect_file(file_path: Path, parent: pytest.Collector):
    if file_path.suffix in {".ipynb"}:
        return NotebookMetadata.from_parent(parent, path=file_path)


def pytest_collection_modifyitems(items: List[pytest.Item]):
    items.sort(key=lambda item: -1 if isinstance(item, KernelSpec) else 0)

    has_oli_api_key: bool = "OLI_API_KEY" in os.environ

    def needs_oli_credentials(item: nbmake.pytest_items.NotebookItem) -> bool:
        # TODO this might have false positives, so we should come up with something better
        # e.g. notebook-level tag
        return "oli" in str(item.filename).lower()

    runtime_tests = (
        item for item in items if isinstance(item, nbmake.pytest_items.NotebookItem)
    )

    for item in runtime_tests:
        if not needs_oli_credentials(item):
            continue
        if not has_oli_api_key:
            item.add_marker(
                pytest.mark.xfail(
                    reason="Notebook needs OLI API credentials, which could not be found in the environment",
                    run=False,
                )
            )
