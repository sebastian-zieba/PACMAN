from pathlib import Path
from typing import Optional


class Logedit:
    """Handles writing output both to terminal and log file."""

    def __init__(self, logname: Path, read: Optional[Path] = None):
        """
        Creates a new log file.

        Parameters
        ----------
        logname : pathlib.Path
            Path of the new log file.
        read : pathlib.Path, optional
            Previous log file whose contents should be copied.
        """

        self.logname = logname

        # Read previous log if requested
        content = []
        if read is not None and read.exists():
            try:
                with read.open("r", encoding="utf-8") as f:
                    content = f.readlines()
            except Exception:
                pass

        # Open new log file
        self.log = self.logname.open("w", encoding="utf-8")

        # Copy previous content
        if content:
            self.log.writelines(content)

    def writelog(self, message, mute=False):
        """
        Print message to terminal and append to log file.
        """

        # Terminal output
        if not mute:
            print(message)

        # Log output
        print(message, file=self.log, flush=True)

    def closelog(self):
        """Close the log file."""
        self.log.close()

    def writeclose(self, message, mute=False):
        """
        Write final message and close log.
        """
        self.writelog(message, mute=mute)
        self.closelog()