import tempfile
import zipfile
import shutil
import os
from typing import Tuple


def update_xlsx_datetime(xlsx_path: str, datetime: Tuple[int, int, int, int,
                                                         int, int]) -> None:
    """Update xlsx file date time, in order to reproduce identically.

    Args:
        xlsx_path: path to the xlsx file
        datetime: a 6-tuple of date time integers
            (year, month, day, hour, minute, second)

    Returns:
        None

    Notes:
        Adapted from the following links:

        - https://stackoverflow.com/a/4653863/4638182
        - https://stackoverflow.com/a/57777044/4638182
    """
    temp_dir_path = tempfile.mkdtemp()
    try:
        temp_xlsx_path = os.path.join(temp_dir_path, "update_datetime.xlsx")
        # For each ZipInfo in the xlsx/zip file infolist()
        # - Read the data of the ZipInfo
        # - Change the date_time of ZipInfo
        # - Write the data and updated ZipInfo to a temp xlsx file
        with zipfile.ZipFile(xlsx_path, "r") as xlsx_read:
            with zipfile.ZipFile(temp_xlsx_path, "w") as xlsx_write:
                for zinfo in xlsx_read.infolist():
                    zdata = xlsx_read.read(zinfo)
                    zinfo.date_time = datetime
                    xlsx_write.writestr(zinfo, zdata)
        shutil.move(temp_xlsx_path, xlsx_path)
    finally:
        shutil.rmtree(temp_dir_path)
