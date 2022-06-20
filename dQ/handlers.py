import logging
import os
import shutil
from abc import ABC
from urllib.parse import unquote
from uuid import uuid4

from redis import Redis
from rq import Queue
from rq.job import Job
from tornado import gen, iostream

import settings
from tornado.web import RequestHandler, stream_request_body

from dQ.operation import Diann

r = Redis(os.getenv("REDIS_HOST"))
q = Queue(connection=r)

class BaseHandler(RequestHandler, ABC):
    def set_default_headers(self):
        self.set_header("access-control-allow-origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', 'GET, PUT, DELETE, OPTIONS')
        self.set_header("Access-Control-Allow-Headers", "access-control-allow-origin,authorization,content-type,unique-id,filename")

    def options(self, **kwargs):
        self.set_status(204)
        self.finish()


class MainHandler(BaseHandler, ABC):
    def get(self):
        self.write(str(uuid4()))


@stream_request_body
class UploadHandler(BaseHandler):
    def initialize(self):
        self.byte_read = 0

    @gen.coroutine
    def data_received(self, chunk):
        if "open_file" not in self.__dict__:
            uniqueID = self.request.headers.get("Unique-ID")
            self.filename = self.request.headers.get("Filename")
            self.folder_path = os.path.join(settings.location, uniqueID)
            os.makedirs(self.folder_path, exist_ok=True)
            os.makedirs(os.path.join(self.folder_path, "temp"), exist_ok=True)
            os.makedirs(os.path.join(self.folder_path, "data"), exist_ok=True)
            self.path = os.path.join(self.folder_path, "temp", self.filename)
            self.open_file = open(self.path, "wb")
        self.open_file.write(chunk)
        self.byte_read += len(chunk)

    def put(self):
        #mtype = self.request.headers.get("Content-Type")
        #logging.info('PUT "%s" "%s" %d bytes', filename, mtype, self.bytes_read)
        self.open_file.close()
        with open(self.path, "rt") as tempfile, \
                open(os.path.join(self.folder_path, "data", self.filename), "wt", newline="") as datafile:
            for line in tempfile:
                line = line.strip()
                if line.startswith("------WebKitFormBoundary"):
                    continue
                elif line.startswith("Content-"):
                    continue
                elif line == "":
                    continue
                else:
                    datafile.write(line + "\n")
        self.write("OK")


def run_diann(folder_path, uniqueID):
    with Diann(os.path.join(folder_path, "data"), os.path.join(folder_path, "DIANN")) as diann:
        print("diann")
    shutil.make_archive(os.path.join(folder_path, uniqueID), "zip", os.path.join(folder_path, "DIANN"))


class DiannHandler(BaseHandler, ABC):
    @gen.coroutine
    def get(self, uniqueID):
        folder_path = os.path.join(settings.location, uniqueID)
        result = q.enqueue(run_diann, job_timeout="2h", args=(folder_path, uniqueID), kwargs={"job_id": uniqueID})
        self.write(result.id)


class ZipHandler(BaseHandler, ABC):
    @gen.coroutine
    def get(self, uniqueID):
        folder_path = os.path.join(settings.location, uniqueID)
        filename = os.path.join(folder_path, uniqueID + ".zip")
        chunk_size = 1024 * 1024 * 1
        with open(filename, 'rb') as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                try:
                    self.write(chunk)
                    yield self.flush()
                except iostream.StreamClosedError:
                    break
                finally:
                    del chunk
                    yield gen.sleep(0.000000001)


class CheckStatusHandler(BaseHandler, ABC):
    @gen.coroutine
    def get(self, uniqueID):
        job = Job.fetch(uniqueID, connection=r)
        print(job.get_status())
        folder_path = os.path.join(settings.location, uniqueID)
        file_path = os.path.join(folder_path, "DIANN", "progress.txt")
        if os.path.exists(file_path):
            with open(file_path, "r") as infile:
                self.write(infile.read())
        else:
            self.set_status(404)
            self.write("Not exists")