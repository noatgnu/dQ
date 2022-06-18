import mimetypes
import os
from functools import partial
from urllib.parse import quote

from tornado import gen
from tornado.httpclient import AsyncHTTPClient
from tornado.testing import AsyncTestCase, gen_test

filenames = ["Reports.log.txt", "Reports.pg_matrix.tsv", "Reports.pr_matrix.tsv", "Reports.stats.tsv"]

@gen.coroutine
def raw_producer(filename, write):
    with open(filename, "rb") as f:
        while True:
            # 16K at a time.
            chunk = f.read(16 * 1024)
            if not chunk:
                # Complete.
                break

            yield write(chunk)


class TestUploader(AsyncTestCase):
    @gen_test
    def test_http_upload(self):
        client = AsyncHTTPClient()
        for filename in filenames:
            mtype = mimetypes.guess_type(filename)[0] or "application/octet-stream"
            headers = {"Content-Type": mtype}
            producer = partial(raw_producer, filename)
            response = yield client.fetch(
                "http://localhost:8000/api/upload/{}".format(filename),
                method="PUT",
                headers=headers,
                body_producer=producer,
            )
            print(response.body)
