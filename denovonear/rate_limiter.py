
import time
import logging
import asyncio

import aiohttp

from denovonear.rate_limiter_retries import ensembl_retry as retry
from denovonear.rate_limiter_retries import ensembl_cache as cache

class RateLimiter:
    ''' class to asynchronously perform http get requests

    This respects the rate limits imposed by the server. Error handling and
    retrying are handled by the 'retry' decorator.
    '''
    MAX_TOKENS = 10
    def __init__(self, per_second=10):
        ''' initialize the class object

        Args:
            per_second: number of queries allowed per second
        '''
        self.tokens = self.MAX_TOKENS
        self.RATE = per_second
        self.updated_at = time.monotonic()
    async def __aenter__(self):
        self.client = aiohttp.ClientSession(trust_env=True)
        return self
    async def __aexit__(self, *err):
        await self.client.close()
        self.client = None

    @cache()
    @retry(retries=9)
    async def get(self, url, headers=None):
        ''' perform asynchronous http get

        Args:
            url: url to get
            headers: http headers to pass in with the get query
        '''
        if not headers:
            headers = {'content-type': 'application/json'}
        await self.wait_for_token()
        async with self.client.get(url, headers=headers) as resp:
            logging.info(f'{url}\t{resp.status}')
            resp.raise_for_status()
            return await resp.read()

    async def wait_for_token(self):
        ''' pause until tokens are refilled
        '''
        while self.tokens < 1:
            self.add_new_tokens()
            await asyncio.sleep(1 / self.RATE)
        self.tokens -= 1

    def add_new_tokens(self):
        ''' add new tokens if sufficient time has passed
        '''
        now = time.monotonic()
        time_since_update = now - self.updated_at
        new_tokens = time_since_update * self.RATE
        if self.tokens + new_tokens >= 1:
            self.tokens = min(self.tokens + new_tokens, self.MAX_TOKENS)
            self.updated_at = now
