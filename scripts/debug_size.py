from metainformant.core.io.download_robust import get_remote_file_size
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14740521/SRR14740521"
print(f"Testing remote size for: {url}")
size = get_remote_file_size(url)
print(f"Size: {size} bytes")

import subprocess
cmd = ["curl", "-s", "-I", "-L", url]
print(f"Running debug curl: {' '.join(cmd)}")
res = subprocess.run(cmd, capture_output=True, text=True)
print("STDOUT:", res.stdout)
print("STDERR:", res.stderr)
