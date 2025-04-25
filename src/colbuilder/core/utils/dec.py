"""
Decorator Utilities for ColBuilder

This module provides utility decorators to enhance functionality and simplify common patterns
in the ColBuilder pipeline. It includes a `timeit` decorator for measuring and logging the
execution time of synchronous and asynchronous functions.

Key Features:
--------------
1. **Execution Time Measurement**:
   - `timeit`: Logs the execution time of a function for performance monitoring.
   - Supports both synchronous and asynchronous functions.

2. **Integration with Logging**:
   - Uses the ColBuilder logging system to log execution times at the DEBUG level.

Usage:
------
This module is designed to be used throughout the ColBuilder pipeline to monitor function
performance and identify bottlenecks.

Example:
--------
```python
from colbuilder.core.utils.dec import timeit

@timeit
def example_function():
    # Simulate some work
    sleep(2)
    return "Done"

@timeit
async def example_async_function():
    # Simulate asynchronous work
    await asyncio.sleep(2)
    return "Async Done"

# Call the functions
example_function()
await example_async_function()
```
"""

# Copyright (c) 2024, Colbuilder Development Team
# Distributed under the terms of the Apache License 2.0

import time
from functools import wraps
from typing import Callable, Any
import asyncio

from colbuilder.core.utils.logger import setup_logger

LOG = setup_logger(__name__)


def timeit(func: Callable[..., Any]) -> Callable[..., Any]:
    """Decorator to measure and log the execution time of a function."""

    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        start_time = time.time()
        if asyncio.iscoroutinefunction(func):

            async def async_wrapper():
                result = await func(*args, **kwargs)
                end_time = time.time()
                LOG.debug(
                    f"{func.__name__} executed in {end_time - start_time:.2f} seconds"
                )
                return result

            return async_wrapper()
        else:
            result = func(*args, **kwargs)
            end_time = time.time()
            LOG.debug(
                f"{func.__name__} executed in {end_time - start_time:.2f} seconds"
            )
            return result

    return wrapper
