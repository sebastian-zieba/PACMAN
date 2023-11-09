import re
from typing import List
from pathlib import Path


def tryint(s):
    try:
        return int(s)
    except Exception:
        return s


def alphanum_key(string: str) -> str:
    """Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]."""
    return [tryint(c) for c in re.split('([0-9]+)', string)]


def sort_nicely(list1: List[Path]):
    """Sort the given list in the way that humans expect."""
    list1.sort(key=lambda x: alphanum_key(x.name))
    return list1
