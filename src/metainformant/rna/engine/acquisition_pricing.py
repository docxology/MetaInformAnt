"""Read-only regional on-demand Linux/gp3 quotes for acquisition scenarios."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime

from metainformant.rna.engine.aws_resources import WorkerPrices, catalog_unit_price


@dataclass(frozen=True, slots=True)
class AcquisitionPriceQuote:
    region: str
    instance_type: str
    disk_gib: int
    observed_at: str
    compute_hourly_usd: float
    gp3_gib_month_usd: float
    hourly_bound_usd: float
    source: str
    disk_throughput_mibps: int = 125
    gp3_mibps_month_usd: float = 0.0


def quote_aws_worker(
    *,
    region: str,
    instance_type: str,
    disk_gib: int,
    profile: str | None = None,
    hourly_floor: float = 0.55,
    disk_throughput_mibps: int = 125,
) -> AcquisitionPriceQuote:
    """Read the catalog; include EBS, one IPv4 address and the existing margin.

    This is on-demand Linux pricing without licenses, CPU-surplus credits,
    nonbaseline gp3 IOPS, or requests/network/object storage. Those
    costs must be reserved separately. The generic launcher uses standard
    CPU credits on burstable instances and rejects licensed AMIs.
    """
    import boto3
    from botocore.config import Config

    pricing = boto3.Session(profile_name=profile, region_name=region).client(
        "pricing",
        region_name="us-east-1",
        config=Config(connect_timeout=10, read_timeout=30, retries={"mode": "standard", "max_attempts": 4}),
    )
    dimensions = {
        "instanceType": instance_type,
        "regionCode": region,
        "operatingSystem": "Linux",
        "tenancy": "Shared",
        "preInstalledSw": "NA",
        "capacitystatus": "Used",
    }
    compute = [
        product
        for page in pricing.get_paginator("get_products").paginate(
            ServiceCode="AmazonEC2",
            Filters=[{"Type": "TERM_MATCH", "Field": k, "Value": v} for k, v in dimensions.items()],
            MaxResults=100,
        )
        for product in page["PriceList"]
    ]
    dimensions = {"volumeApiName": "gp3", "regionCode": region, "productFamily": "Storage"}
    storage = [
        product
        for page in pricing.get_paginator("get_products").paginate(
            ServiceCode="AmazonEC2",
            Filters=[{"Type": "TERM_MATCH", "Field": k, "Value": v} for k, v in dimensions.items()],
            MaxResults=100,
        )
        for product in page["PriceList"]
    ]
    throughput_price = 0.0
    if disk_throughput_mibps > 125:
        dimensions = {"volumeApiName": "gp3", "regionCode": region, "productFamily": "Provisioned Throughput"}
        throughput = [
            product
            for page in pricing.get_paginator("get_products").paginate(
                ServiceCode="AmazonEC2",
                Filters=[{"Type": "TERM_MATCH", "Field": k, "Value": v} for k, v in dimensions.items()],
                MaxResults=100,
            )
            for product in page["PriceList"]
        ]
        throughput_price = catalog_unit_price(throughput, "GiBps-mo") / 1024
    prices = WorkerPrices(
        catalog_unit_price(compute, "Hrs"),
        catalog_unit_price(storage, "GB-Mo"),
        hourly_floor,
        gp3_mibps_month=throughput_price,
    )
    return AcquisitionPriceQuote(
        region,
        instance_type,
        disk_gib,
        datetime.now(UTC).isoformat(),
        prices.compute_hourly,
        prices.gp3_gib_month,
        prices.hourly_bound(disk_gib, throughput_mibps=disk_throughput_mibps),
        (
            "AWS GetProducts on-demand Linux and selected gp3 throughput at baseline IOPS; "
            "conservative 28-day storage month"
        ),
        disk_throughput_mibps,
        throughput_price,
    )
