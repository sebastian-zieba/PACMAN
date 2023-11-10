import re
from typing import List, Union
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


def sort_nicely(input_list: List[Union[str, Path]]) -> List[Union[str, Path]]:
    """Sort the given list in the way that humans expect."""
    return sorted(input_list, key=lambda x: alphanum_key(
        x.name if isinstance(x, Path) else x))
