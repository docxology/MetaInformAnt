"""Explicit platform and capacity policy for the default acquisition bootstrap."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class WorkerImage:
    state: str
    architecture: str
    platform_details: str
    has_product_codes: bool


def validate_worker_image(image: WorkerImage, *, custom_template: bool = False) -> None:
    """Reject platforms whose software or license costs the model cannot cover."""
    if image.state != "available" or image.platform_details != "Linux/UNIX" or image.has_product_codes:
        raise ValueError("generic acquisition requires an available, unlicensed Linux AMI")
    if not custom_template and image.architecture != "x86_64":
        raise ValueError(
            "default acquisition bootstrap requires x86_64; other architectures need a custom startup template"
        )
