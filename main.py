import asyncio
import os
import sys
import tornado.httpserver
from tornado.ioloop import IOLoop
from tornado.web import Application
from tornado.options import define, options

from dQ.handlers import MainHandler, UploadHandler, DiannHandler, ZipHandler, CheckStatusHandler

define("port", default=8000, help="Port number")

if sys.platform.startswith("win32"):
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

routes = [
    (r"/", MainHandler),
    (r"/api/upload/", UploadHandler),
    (r"/api/diann/(.*)/fasta/(.*)/gaf/(.*)/obo/(.*)/", DiannHandler),
    (r"/api/download/(.*)/", ZipHandler),
    (r"/api/status/(.*)/", CheckStatusHandler),
]

settings = {
    "debug": False,
    "autoreload": False,
    "autoescape": True,
    "x-header": True
}

if __name__ == '__main__':
    tornado.options.parse_command_line()
    app = Application(routes, **settings)
    app.listen(options.port)
    IOLoop.current().start()
