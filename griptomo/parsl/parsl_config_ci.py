import parsl
from parsl.config import Config
from parsl.executors.threads import ThreadPoolExecutor


def make_ci_config():
    return Config(
        executors=[ThreadPoolExecutor(max_threads=4, label="threads")],
        strategy="none",
    )
