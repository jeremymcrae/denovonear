
from pathlib import Path
import asyncio
import random
import functools

import aiohttp

from denovonear.ensembl_cache import EnsemblCache

def ensembl_cache(folder=None):
    ''' store/retrive repeated ensembl requests from a persistent sqlite cache
    '''
    if folder is None:
        folder = Path.home() / '.cache' / 'ensembl'
    cache = EnsemblCache(str(folder))
    def decorator(func):
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            url = args[1]
            if 'rest.ensembl.org' in url:
                cached = cache.get_cached_data(url)
                if cached is not None:
                    return cached
                data = await func(*args, **kwargs)
                cache.cache_url_data(url, data)
                return data
            else:
                return await func(*args, **kwargs)
        return wrapper
    return decorator

def ensembl_retry(retries=5):
    ''' perform all the error handling for the request
    
    retries up to N times under certain error conditions, and increases waiting
    time between requests, unless we've hit rate limits, when it uses the stated
    retry time.
    '''
    def decorator(func):
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            result = None
            last_exception = None
            for i in range(retries):
                try:
                    return await func(*args, **kwargs)
                except (aiohttp.ServerDisconnectedError, aiohttp.ClientOSError,
                        asyncio.TimeoutError) as err:
                    last_exception = err
                    delay = 0
                except aiohttp.ClientResponseError as err:
                    last_exception = err
                    # 500, 503, 504 are server down issues. 429 exceeds rate
                    # limits. 400 is server memory issue. Raises other errors.
                    if err.status not in [500, 503, 504, 429, 400]:
                        raise err
                    if err.status == 400:
                        if 'Cannot allocate memory' not in err.message:
                            raise err
                    delay = random.uniform(0, 2 ** (i + 2))
                    if err.status == 429:
                        delay = float(dict(err.headers)['Retry-After'])
                await asyncio.sleep(delay)
            if last_exception:
                raise last_exception
            return result
        return wrapper
    return decorator
