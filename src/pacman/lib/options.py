import sys

OPTIONS = {}

if sys.platform == "win32":
    OPTIONS["encoding"] = "utf-8-sig"
else:
    OPTIONS["encoding"] = "utf-8"
